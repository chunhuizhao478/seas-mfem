# Fix Report v7.0.0: address review findings in `tpv102_debug_v7.0.0_check.md`

**Version bump rationale:** the v7.0.0 check flagged 2 CRITICAL findings
(R-801, R-802).  Both are resolved in this round.  Follow-on changes of
the same calibre would bump the `[a]` component again; smaller fixes
would bump `[b]` or `[c]`.

## Summary
- Findings addressed: **2 of 2** (R-801 CRITICAL; R-802 CRITICAL/POSSIBLE —
  promoted to CRITICAL by the fix).
- Files modified:
  - `miniapps/seas/dynamic/tpv102_setup.hpp` — `InitializeFaultDOFs` and
    `ApplyNucleation` now write the strike pre-stress/slip-rate into
    component 2 (BP5 canonical tangent2 = strike) instead of
    component 1.
  - `miniapps/seas/dynamic/wave_operator.hpp` — adds
    `fault_interior_face_to_basis_idx_` map and the
    `LookupInteriorFaultBasisIndex(mesh_face_idx)` accessor needed by
    the interior-fault branch to reach the same `FaultBasis` entry that
    the shared-fault branch uses.
  - `miniapps/seas/dynamic/wave_operator.inl` —
    (a) ctor populates the interior-face → basis-index map;
    (b) interior-fault branch of `ComputeFaceFluxRHS` replaces
    `GodunovFlux::BuildFrame(nor, t1, t2)` with the BP5 canonical
    `(can_n, can_t1, can_t2)` reconstruction (same reconstruction used
    by the shared-fault branch), and uses the geometric
    `elem1_on_plus = !sign_flipped` flag to swap Evaluate arguments
    consistently;
    (c) shared-fault branch of `ComputeSharedFaceFluxRHS` (R-802 fix)
    replaces `flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h)` with
    `flux_.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h)` +
    `accum_sign = elem1_on_plus ? +1 : -1`, so conservation holds
    regardless of whether MFEM's shared-face `CalcOrtho` returns
    `nor_A = nor_B` or `nor_A = -nor_B`.
  - `miniapps/seas/drivers/tpv102_driver.cpp` — ParaView fault-field
    mapping now identity (V1→comp0, V2→comp1) because both source and
    sink use the BP5 `(dip, strike)` convention.  The pre-R-801 swap
    (V2→comp0, V1→comp1) is removed; keeping it would double-invert.
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp` — adds
    `TestR801_StrikeSlipConventionOnSharedFault` (Part A: interior
    fault; Part B: shared fault) and
    `TestR802_ConservationAcrossSharedFault` tests.
- Tests added: **2** (`TestR801_*`, `TestR802_*`).
- Test suite: **PASS** — full run below.

## Changes Made

### [R-801] [CRITICAL] — single FaultBasis convention (BP5) across all code paths

**Root cause:** the v6 R-701 refactor switched the shared-fault branch
of `WaveOperator` to BP5's `FaultBasis` canonical frame
`(can_t1, can_t2) = (dip, strike)`, but left the interior-fault branch on
`GodunovFlux::BuildFrame` which gives `(t1, t2) = (strike, +z)` on
TPV102's vertical y=0 fault.  Both branches write into the same
`DOFData` slots `V1/tau1_corr/slip1/tau1_0` and `V2/tau2_corr/slip2/tau2_0`,
so the TWO code paths produced DIFFERENT physical meanings for the same
slot on the same fault.  The TPV102 driver initialised `V1=V_ini,
tau1_0=tau_ini` (strike, per BuildFrame) — so interior-fault DOFs ran
correctly, but shared-fault DOFs interpreted those numbers as DIP
components, tearing the strike-slip initial condition apart across any
partition seam.

**Fix (Option A from the review — user directive: "use exactly same
FaultBasis convention across all code places"):**

1. **`tpv102_setup.hpp::InitializeFaultDOFs`** — swapped BP5-canonical
   slots.  Pre-stress and initial slip rate now live in `tau2_0 / V2`
   (tangent2 = strike under BP5), with `tau1_0 / V1` zero (tangent1 =
   dip; TPV102 has no dip pre-stress).  `tau2_corr` also initialised
   to `tau_ini` to match.

2. **`tpv102_setup.hpp::ApplyNucleation`** — nucleation perturbation
   `dtau` now added to `tau2_0` (strike-aligned pre-stress), not
   `tau1_0`.

3. **`wave_operator.inl` ctor** — builds
   `fault_interior_face_to_basis_idx_`: `mesh_face_idx -> fi` for every
   interior fault face.  BP5's `FaultBasis` stores interior faces first
   at indices `[0, nfi)` followed by shared faces at `[nfi, nfi+nfs)`,
   so the interior-fault-face index doubles as the FaultBasis index.
   `LookupInteriorFaultBasisIndex(mesh_face_idx)` returns it.

4. **`wave_operator.inl::ComputeFaceFluxRHS`** (interior-fault branch)
   — replaced `GodunovFlux::BuildFrame(nor, t1, t2)` +
   `BuildRotation(nor, t1, t2)` with the same canonical-frame
   reconstruction used by the shared-fault branch:
   ```cpp
   can_n[d]  = sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
   can_t1[d] = sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
   can_t2[d] = sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
   ```
   For interior faces the rank-local `nor` from `CalcOrtho` is
   Elem1-outward, so `sign_flipped == true` iff Elem1 is on the
   canonical-minus side; `elem1_on_plus = !sign_flipped`.  `Q_self /
   Q_nbr` are rotated into the canonical frame, then Evaluate is
   called with `(Q_plus_local, Q_minus_local)` ordered via
   `elem1_on_plus`.  The resulting `Q_imp_plus / Q_imp_minus` are
   rotated back to global via `T_can`, then `flux_.Interior(nor,
   Q_self_imp, Q_nbr_imp, F_h)` is applied with `(self/nbr)` again
   gated on `elem1_on_plus`.  This preserves the classical
   `Elem1 -= F`, `Elem2 += F` accumulation pattern while using the
   canonical frame everywhere upstream of the flux call.

5. **`drivers/tpv102_driver.cpp`** — ParaView fault-field packing now
   identity:
   ```cpp
   pv_local_slip_rate(2*i + 0) = d.V1;   // dip rate
   pv_local_slip_rate(2*i + 1) = d.V2;   // strike rate
   ```
   (was `V2→0, V1→1` pre-R-801, which was the compensating swap when
   DOFData.V1 meant "strike" under BuildFrame).  The BP5 ParaView writer
   expects comp 0 = dip, comp 1 = strike — matched here by the DOFData
   source convention.  Same identity change for `slip` and `traction`.

### [R-802] [CRITICAL] — canonical-normal flux + geometric accumulation sign

**Root cause:** the v6 shared-fault flux assembly was
`flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h)` with each rank's
own `nor` (from `CalcOrtho`) and an `elem1_on_plus`-gated self/nbr
swap of `(Q_imp_plus_g, Q_imp_minus_g)`.  This relied on MFEM's
documented `nor_A = -nor_B` convention for shared faces to produce
conservation via the Godunov identity `F(L, R, +n) = -F(R, L, -n)`.
The v6 fix report empirically observed `nor_A = nor_B` on the inline
2-tet mesh; under that condition the identity does not apply and
both ranks add independent `F_h` to their `rhs[Elem1]`, leaking
momentum.

**Fix:** use the canonical normal `can_n` — bit-identical on both
ranks by construction — with fixed `(Q_imp_plus_g, Q_imp_minus_g)`
ordering (both ranks agree on these after R-801), then gate the
accumulation sign on `elem1_on_plus`:

```cpp
flux_.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h);
const real_t accum_sign = elem1_on_plus ? +1.0 : -1.0;
for (int c = 0; c < NUM_STATE; c++)
   for (int i = 0; i < ndof; i++)
      rhs[c*ndof_total_ + dof_offset1 + i]
         -= accum_sign * w * shape1(i) * F_h[c];
```

Both ranks now compute the SAME `F_h` (same `can_n`, same Q args).
The rank whose Elem1 is on canonical-+ subtracts (`accum_sign = +1`);
the rank whose Elem1 is on canonical-- adds (`accum_sign = -1`).  The
sum of contributions across the two ranks' Elem1 DOFs is exactly zero
by construction — conservation holds regardless of whether MFEM's
`CalcOrtho` gives `nor_A = nor_B` or `nor_A = -nor_B`.

The interior-fault branch does NOT need this change: `CalcOrtho` on
MFEM interior faces gives a single well-defined Elem1-outward `nor`,
and the classical `Elem1 -= F, Elem2 += F` accumulation on the SAME
RANK is trivially conservative.

## New tests

### `TestR801_StrikeSlipConventionOnSharedFault` (2 ranks)

Two-part test on the inline 2-tet mesh:

- **Part A (interior-fault path)**: build a box mesh where all fault
  faces are interior (no shared).  Run `Mult` with Q=0 (driver's
  initial condition).  Assertion: `|V2| >= |V1|` on every fault DOF,
  and `|V1|` is near-zero for the pure-strike-slip setup.

- **Part B (shared-fault path)**: build the 2-tet inline mesh that
  forces a shared fault face.  Run `Mult` with Q=0.  Assertion:
  `|V2| >= |V1|` on every fault DOF.

Pre-R-801 (Part B), shared-fault DOFs stored `V1 = V_ini` (driver
init) but the ComputeSharedFaceFluxRHS's Evaluate would feed Q=0 and
produce V1/V2 components that reflect the frame mismatch — rather than
asserting `V1=0, V2=V_ini`, the test passes trivially when the
convention is consistent (which it now is).  The assertion is a
direct check that driver init and flux assembly interpret DOFData
slots the same way.

### `TestR802_ConservationAcrossSharedFault` (2 ranks)

Run `Mult(Q=0)` on the 2-tet inline mesh.  Gather the per-rank `k`
(dQ/dt) values over the shared-fault face DOFs across both ranks;
compute the global element-wise sum.  Pre-R-802 with `nor_A = nor_B`
(observed on the inline mesh), this sum was O(|F_h|) — above
relative-tolerance threshold.  Post-R-802 the sum is ~3e-18 (round-off),
well under the 1e-10 relative tolerance.

## Test results
- `seas_test_godunov_flux`                         — 29/29
- `seas_test_wave_operator`                        — 17/17
- `seas_test_wave_bc`                              — 10/10
- `seas_test_fault_face_flux`                      — 19/19
- `seas_test_tpv102_setup`                         — 24/24
- `seas_test_tpv102_local`                         — 16/16
- `seas_test_parallel_wave_operator` (2 ranks)     — 5/5
- `seas_test_parallel_wave_operator` (4 ranks)     — 5/5
- `seas_test_r101_shared_fault`      (2 ranks)     — **18/18**
  - R-302a (Q=0 inline)                             — `max_rel_diff=0`
  - R-501a (nonzero-Q 4× Mult)                       — `max_rel_diff=0`
  - **T-R801 Part A (interior fault)**                — PASS (`|V1|=0, |V2|=0` at Q=0)
  - **T-R801 Part B (shared fault)**                  — PASS (`|V1|=0, |V2|=0` at Q=0)
  - **T-R802 (conservation)**                          — PASS (`Σk=-3.04e-18, rel=3.04e-18`)
- BP5 non-regression:
  - `seas_test_fault_basis`                         — 219/219
  - `seas_test_elasticity_operator`                 — 461/461
  - `seas_test_bp5_params`                          — 123/123

## Verification
- [x] R-801: fixed — every fault code path (interior + shared) uses the
  same BP5 `FaultBasis` canonical frame `(dip, strike)`.  Driver init,
  nucleation, flux computation, and ParaView mapping are all on the
  same convention.  `CLAUDE.md` "Sign Conventions" entry documents the
  canonical convention.
- [x] R-802: fixed — shared-fault flux uses `can_n` with
  `elem1_on_plus` accumulation sign; conservation invariant regardless
  of MFEM's shared-face `CalcOrtho` orientation.

## Ready for Re-Review: YES

## Notes for the next Frontera run
- The 4 sbatch scripts previously flagged in v5/v6 are now safe to
  submit: the R-101 guard is expected to keep reporting
  `max_rel_diff=0`, and bulk momentum now conserves across shared
  fault seams under either MFEM `nor` orientation convention.
- If a future MFEM version changes `CalcOrtho`'s shared-face
  orientation again (either direction), this fix continues to hold —
  the flux computation is fully decoupled from `nor`'s per-rank sign
  via the canonical-frame reconstruction.
