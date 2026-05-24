# Reference comparison: tension / no-opening / normal-channel handling — ours vs SeisSol vs DRDG3D — 2026-05-23

## Scope
Compare how three dynamic-rupture codes treat the fault NORMAL channel under tension, to settle
whether the SAFS tensile-σ_n runaway is a tension-handling defect in our code. **No code modified.**
Reference trees: `/Users/chunhuizhao/projects/SeisSol` (C++), `/Users/chunhuizhao/projects/drdg3d`
(Fortran 90). Ours: `dynamic/fault_face_flux.cpp` (`ComputeTrialTraction`, `BuildImposedState`),
`dynamic/tpv205_friction.hpp` (`SolveLSW_TPV205`), `dynamic/tpv205_substep_iterator.cpp`.

**Verification note:** SeisSol lines and the DRDG3D slip-metric + strength lines were read and
verified directly. DRDG3D's flux internals (`fstar`) and its rate-state slip-rate limiter are from
the sub-agent read and are marked [unverified-here] where I did not re-open them.

## Bottom line
On the two rules that matter for the runaway — **(1) clamp the strength under tension via
`max(σ_n,0)` and (2) weld the normal channel (no opening): our code matches SeisSol term-for-term
and matches DRDG3D in effect.** None of the three codes bounds the transmitted tensile normal
traction, and neither our nor SeisSol's slip-weakening path caps the slip rate. So the runaway is
**not** a tension-handling deviation — it is the bulk producing spurious *sustained* tension, which
the reference codes avoid by resolution (and SeisSol additionally by *resampling* the slip rate, a
measure we lack). Two of our items differ from SeisSol but are defensible (DRDG3D agrees with us on
one of them).

---

## Findings (comparison points, ranked by relevance to the runaway/speckle)

### [C-001] CONFIRMED MATCH — no-opening / welded normal channel: our `BuildImposedState` is structurally identical to SeisSol's imposed-state, and to DRDG3D in effect

**Relevance:** decisive — this is the rule the runaway was suspected to violate.

SeisSol (`FrictionLaws/FrictionSolverCommon.h:347-361`, verified):
```cpp
imposedStateM[N][i] += weight * normalStress;                 // transmit Godunov normal stress
imposedStateM[U][i] += weight * (qIMinus[U] - invZpNeig*(normalStress - qIMinus[N]));  // normal vel
...
imposedStateP[N][i] += weight * normalStress;
imposedStateP[U][i] += weight * (qIPlus[U]  + invZp    *(normalStress - qIPlus[N]));
```
Ours (`fault_face_flux.cpp:278,286,291,295`):
```cpp
Q_imp_minus[VX] = Q_minus[VX] - invZp_m*(s.sigma_n_corr - Q_minus[SXX]);
Q_imp_plus[VX]  = Q_plus[VX]  + invZp_p*(s.sigma_n_corr - Q_plus[SXX]);
Q_imp_*[SXX]    = s.sigma_n_corr;     // = sigma_n_trial (Godunov normal stress, transmitted)
```
These are the **same formula** (`N↔SXX`, `U↔VX`, `normalStress↔σ_n_corr`). In both, friction sets
only the two SHEAR tractions (`traction1/2` ↔ `tau1/2_corr`); the normal stress is the Godunov
trial value passed through **unchanged**. The imposed normal-velocity jump is therefore identically
zero in both (the algebra: `v_n⁺_imp − v_n⁻_imp = 0`), i.e. **welded / no-opening** — exactly TPV14
rule 2. DRDG3D welds likewise (`mod_wave.F90:680-688`, normal flux from the Godunov state, only
shear modified) [unverified-here].

**Conclusion:** our no-opening handling is correct and matches both references. Neither reference
clamps the transmitted normal stress — both transmit the raw (possibly tensile) Godunov value, as
we do. So the −9 TPa is a welded interface faithfully transmitting a spurious *bulk* tensile stress
in all three codes' formulation; it is not an opening artifact and not a deviation.

---

### [C-002] CONFIRMED MATCH — tension strength clamp `max(σ_n,0)` is identical in all three

**Relevance:** decisive — TPV14 rule 1.

- Ours (`tpv205_friction.hpp:133-134`): `sigma_n_pos = max(sigma_n_total,0); tau_strength = mu_eff*sigma_n_pos;`
- SeisSol (`LinearSlipWeakening.h:158-160`, verified): `strength = -cohesion - mu*std::min(totalNormalStress, 0.0);` — compression is **negative** in SeisSol, so `min(σ,0)` keeps compression and zeroes tension. Same effect as our `max(σ,0)` in our compression-positive convention.
- DRDG3D (`mod_wave.F90:837`, verified): `Tau_str = fault_mu*max(0.0,-Tau_n) + C0;` — same clamp.

**Conclusion:** identical across all three. Our rule-1 handling is correct. (SeisSol and DRDG3D also
carry a cohesion `C0`/`cohesion` term; ours folds prestress differently but the σ_n clamp is the
same.)

---

### [C-003] DIFFERENCE — SeisSol RESAMPLES the slip rate before integrating slip (anti-aliasing of grid-scale oscillations); we do not

**Relevance:** HIGH for the speckle (grid-scale ringing), lower for the σ_n runaway.

SeisSol (`LinearSlipWeakening.h:174-186`, verified):
```cpp
real resampledSlipRate[...];
specialization.resampleSlipRate(resampledSlipRate, this->slipRateMagnitude[ltsFace]);
...
const auto update = resampledSlipRate[pointIndex] * this->deltaT[timeIndex];
this->accumulatedSlipMagnitude[...] += update;     // integrate the RESAMPLED rate
```
with the comment (`:227-228`): *"Resample slip-rate, such that the state increment (slip) lies in
the same polynomial space as the degrees of freedom."* This is the Pelties et al. projection that
suppresses spurious high-frequency (grid-scale) slip-rate oscillations in DG dynamic rupture — the
exact class of artifact as the SAFS "speckle." **We integrate the raw slip rate with no resample**
(`tpv205_substep_iterator.cpp:122-123`).

**Caveat / honesty:** SeisSol's effect is largest at p≥2; SAFS runs at **p=1** (plan B.2.5), where
the SSO deficit is weak, so resampling would help less here than at high order. Still, it is a
genuine reference technique aimed squarely at grid-scale speckle and worth noting as a candidate
mitigation if the speckle survives resolution refinement.

**Action (not a fix to apply now):** if the resolution study (R-007) leaves residual grid-scale
oscillation, a SeisSol-style slip-rate resample/projection is the reference-blessed remedy. Out of
scope for the tension question.

---

### [C-004] DIFFERENCE — no upper cap on slip rate in our SW path nor SeisSol's; DRDG3D has a limiter only in its rate-state path

**Relevance:** MODERATE — bears on the proposed `V≫c_s` sanity bound (debug-review R-002).

- Ours (`tpv205_friction.hpp:137-144`): `V_abs = max(0,(tau_abs-tau_str)/eta_s)` — non-negative, **no upper bound**.
- SeisSol (`LinearSlipWeakening.h:77-79`, verified): `slipRateMagnitude = max(0,(absoluteTraction-strength)*invEtaS)` — non-negative, **no upper bound**.
- DRDG3D (`mod_wave.F90:878-879`, [unverified-here]): rate-state path only — `if (V > Phi) V = 0.5*Phi/eta;` a stress-derived slip-rate limiter. Not present in its slip-weakening path.

**Conclusion:** our lack of an upper V cap in the LSW path **matches SeisSol** (the closest analogue
to our LSW). A `V≫c_s` sanity bound would therefore be a *beyond-SeisSol* numerical safeguard — it
has a partial precedent in DRDG3D's rate-state limiter but is not what either reference does for
slip-weakening. It is defensible as a debugging tripwire (a slip rate ≫ c_s is unphysical), but
should be presented as our addition, not as matching the references.

---

### [C-005] DIFFERENCE — slip-weakening distance: SeisSol uses PATH LENGTH (∫|V|dt, TPV14-strict); we and DRDG3D use NET slip (|∫V dt|)

**Relevance:** LOW for the runaway (monotonic slip ⇒ net ≈ path); relevant only under backwards motion.

- TPV14 spec: `d` = total path length (the screenshot's "backwards motion" paragraph).
- SeisSol (`LinearSlipWeakening.h:185-186 + 234-239`, verified): `accumulatedSlipMagnitude += resampledSlipRate*dt` (the slip-rate **magnitude** integrated ⇒ ∫|V|dt = path length), used as `min(fabs(accumulatedSlip)/Dc, 1)`. **Path length — TPV14-compliant.**
- DRDG3D (`seis3d.F90:319-352`, verified): `slip1 += rk*tslip1; slip = sqrt(slip1²+slip2²)` ⇒ `|∫V dt|` = **NET** slip. (The sub-agent's "monotonic path length" label was wrong; the formula is the net magnitude.)
- Ours (`fault_face_flux.cpp:794`, `tpv205_substep_iterator.cpp:99`): `delta = sqrt(slip1²+slip2²)` = **NET** slip — same as DRDG3D.

**Conclusion:** the references **disagree** here — SeisSol follows TPV14 (path length), DRDG3D
matches us (net). For the runaway DOF (≈monotonic slip) the two coincide, so this is not the cause.
It is a genuine TPV14-compliance gap vs SeisSol that matters only under oscillatory/reversing slip;
switching to path length would *weaken oscillatory cells faster* (more slip), so it is not a speckle
fix and should be treated purely as an optional spec-compliance item.

---

## Summary
- Decisive matches (ours ≡ references): **C-001 no-opening/welded normal channel** (identical to
  SeisSol term-for-term), **C-002 `max(σ_n,0)` strength clamp** (identical in all three). Our
  tension handling is correct and matches the references.
- Differences: **C-003** SeisSol resamples the slip rate (anti-SSO) — we don't (matters most at
  p≥2; SAFS is p=1); **C-004** no upper V cap in our or SeisSol's SW path (a `V≫c_s` bound would be
  our addition, partial precedent in DRDG3D rate-state); **C-005** SW distance is path-length in
  SeisSol (TPV14) vs net in ours/DRDG3D (only matters under backwards motion).
- **Verdict on the question "how should we handle tension / are we doing it wrong":** we are doing
  it the same as SeisSol and DRDG3D. The tensile-σ_n runaway is **not** a tension-handling bug — all
  three codes transmit the raw tensile Godunov normal stress through a welded fault and clamp only
  the strength. The references avoid the runaway by **resolution** (a well-resolved cohesive zone
  rarely sustains tension — TPV14: "not likely to encounter tension") and, in SeisSol, by
  **slip-rate resampling**. This triangulates with the earlier reviews: the lever is upstream
  (process-zone resolution / the dip→normal feedback), optionally plus SeisSol-style slip-rate
  resampling for grid-scale oscillation and a beyond-reference `V` sanity bound for debugging.

## Recommended next steps (no code change in this report)
1. **Treat the tension handling as settled/correct** — do not alter `max(σ_n,0)` or the welded
   normal channel; they match SeisSol exactly (C-001/C-002).
2. **Resolution study (R-007)** remains the primary lever for the runaway/speckle.
3. If grid-scale speckle persists after refinement, **port SeisSol's slip-rate resample** (C-003) —
   the reference-standard anti-SSO measure — as its own task; note the p=1 caveat.
4. The **`V≫c_s` sanity tripwire** (debug-review R-002) is fine as a *debugging* guard but is our
   addition (C-004), not a reference behavior; keep it non-fatal/diagnostic.
5. **C-005 (path-length d)** is an optional TPV14-compliance item vs SeisSol; not a runaway fix.

## Unreviewed / lower-confidence areas
- DRDG3D's `get_flux` normal-channel internals (`fstar(1)/fstar(4)`) and its rate-state slip-rate
  limiter (`mod_wave.F90:680-688, 878-879`) are from the sub-agent read; not re-opened line-by-line
  here. The DRDG3D slip-metric (net) and SW strength (`max(0,-Tau_n)`) WERE verified directly.
- SeisSol's pore-pressure / `fluidPressure` term in `totalNormalStress` (`LinearSlipWeakening.h:156`)
  is present but irrelevant to the dry SAFS comparison; not analyzed.
- Whether SeisSol's `resampleSlipRate` default specialization is the identity (no-op) or the
  projection depends on the friction-law specialization selected at build; the NoSpecialization
  path (`:232`) does project, the copy path (`:266,314`) does not — confirm which TPV-class config
  uses which before porting (C-003).
