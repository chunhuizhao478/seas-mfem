# B0 verdict: per-side bi-material FAULT Riemann verification

Date: 2026-06-06
Phase: Part B / B0 (GATE) of `document/bimaterial_dev/PLAN_unified_bimaterial_volume_and_fault_2026-06-06.md`
Test: `tests/unit/test_bimaterial_fault_riemann.cpp` (`seas_test_bimaterial_fault_riemann`)

## VERDICT: PASS — per-side fault flux is correct for the matrix path.

B2 (guard relaxation) is UNBLOCKED.  The fault-flux math
(`FaultFaceFlux::ComputeTrialTraction` / `BuildImposedState`) already consumes per-side
`Zp_plus/Zp_minus`, `Zs_plus/Zs_minus`, `eta_p/eta_s`, and is the CORRECT bi-material
welded Riemann for UNEQUAL impedance (not merely syntactically per-side).

## Analytic verification (test 10/10, machine precision)

Welded-interface star (derived; matches `ComputeTrialTraction` algebraically, so the
match is sign-convention-independent):

    sigma* = (Z2*sigma_L + Z1*sigma_R + Z1*Z2*(v_L - v_R)) / (Z1 + Z2)
    v*     = (Z1*v_L + Z2*v_R + (sigma_L - sigma_R)) / (Z1 + Z2)   [continuous]

- `sigma_n_trial == sigma*` (P, per-side Zp) and `tau{1,2}_trial == tau{1,2}*`
  (S, per-side Zs) for a general bimaterial state — rel <= 2e-16.
- Right-going incident P-wave (sigma_i = +Z1*v_i, side 2 at rest) reproduces the
  welded R/T EXACTLY: stress reflection `R_sigma = (Z2-Z1)/(Z1+Z2) = 0.3151`,
  velocity transmission `T_v = 2*Z1/(Z1+Z2) = 0.6849` (TPV6 contrast Z1=Zp(far),
  Z2=Zp(near)).
- `BuildImposedState` (locked, V=0): imposed normal velocity is CONTINUOUS across the
  fault (`v_imp_plus == v_imp_minus == T_v*v_i`) and the imposed traction is
  single-valued — i.e. the locked fault reproduces the welded interface.
- Homogeneous reduction (Z1==Z2): `sigma_n_trial == (sigma_L+sigma_R)/2 + Z(v_L-v_R)/2`
  (the standard homogeneous Riemann) — confirms byte-exact reduction for
  TPV31/TPV102/BP5 (where Zp_plus==Zp_minus).

## Conversion-site TABLE (R-005): every guarded variant's imposed-state -> flux site is per-side A

The OPERATOR (Mult / AdvanceADER), not the Evaluate variant, performs the
imposed-state -> bulk-flux conversion.  On the `BimaterialWaveOperator`, ALL three
conversion sites apply PER-SIDE A via the virtual `FluxForElem_` (verified by reading
`wave_operator.inl`):

| conversion site (wave_operator.inl) | code                                              | per-side A? |
|-------------------------------------|---------------------------------------------------|-------------|
| Mult / RK, interior fault face `:3048-3051` | `FluxForElem_(elem_plus).Interior(...)` + `FluxForElem_(elem_minus).Interior(...)` | YES |
| ADER, interior fault face `:4301-4303`      | `FluxForElem_(elem_plus).Interior(...)` + `FluxForElem_(elem_minus).Interior(...)` | YES |
| shared (cross-rank) fault face `:3715`      | `FluxForElem_(e1).Interior(...)` (each rank applies ITS local side's A)            | YES |

`BimaterialWaveOperator` does NOT override `Mult`, so it inherits these per-side
conversions (verified: no `void Mult` in `bimaterial_wave_operator.hpp`).

Guarded variants (`fault_face_flux.cpp`) and the operator path that dispatches them:

| variant (guard line)                      | dispatched by            | conversion | matrix-path target            |
|-------------------------------------------|--------------------------|------------|-------------------------------|
| `Evaluate` (:359)                         | Mult/RK (rate-state)     | :3048 per-side | spatial RS; native TPV102/104 |
| `EvaluateTotal` (:499)                    | Mult/RK (total-Q RS)     | :3048 per-side | TPV102 total-Q (scalar today) |
| `EvaluateADER_LSW` (:778)                 | AdvanceADER (LSW)        | :4301 / :3715 per-side | TPV6 Arm 1 (ADER), TPV31 |
| `EvaluateLSW` (:912)                      | Mult/RK (LSW)            | :3048 per-side | TPV6 Arm 2 (RK), TPV31 mixedflux |
| `EvaluateADER_LSW_ForcedRupture` (:1030)  | AdvanceADER (LSW forced) | :4301 per-side | TPV26/27 forced rupture       |

Because the conversion is per-side at every operator site, and the relaxation flag is
set ONLY by the `BimaterialWaveOperator` (B2), relaxing all five guards via the flag is
safe: any variant the matrix operator dispatches converts per-side.  The scalar
`WaveOperator` never sets the flag, so its single-A conversion remains guarded (aborts
on a bimaterial fault).

## Conclusion
No defect found.  Proceed to B1 (per-side material assignment) and B2 (guard
relaxation, matrix-operator flag).
