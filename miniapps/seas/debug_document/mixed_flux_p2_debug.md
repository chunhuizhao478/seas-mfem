# p=2 Mixed-Flux Adjacent — Catastrophic Blow-up Debug Report

**Date:** 2026-04-29
**Reporter:** debug
**Affected runs:**
- `tpv102/results_mixed_flux_adjacent_p2_O3_normal_job7681805`
- `tpv104/results_mixed_flux_adjacent_p2_O3_normal_job7680797`
- `tpv205/plots_results_mixed_flux_adjacent_p2_O3_normal_job7682076`

---

## Symptom

| Test | Mode | Polynomial / ADER | Result |
|------|------|-------------------|--------|
| TPV102 | Adjacent | p=2 / O3 | Slip-rate oscillations, ringing in σ_n, slip-strike erosion across propagated area |
| TPV104 | Adjacent | p=2 / O3 | Slip rate blows up to **~1e11 m/s** at t≈7.8s |
| TPV205 | Adjacent | p=2 / O3 | Slip rate blows up to **~1e10 m/s** at t≈6.5s |
| TPV102 | Adjacent | **p=1** / O2 | Works |
| TPV102 | None    | **p=2** / O3 | Works |

The user reports the systematic erosion / isolated high values appear **across the propagated area, not just at peaks** — consistent with an instability that grows wherever a wave passes a fault-adjacent (central-flux) face.

---

## Reproducer (built and run)

`miniapps/seas/tests/unit/test_mixed_flux_dispatch_adjacent_p2.cpp` — parallel to the existing `test_mixed_flux_dispatch_adjacent.cpp` (p=1) but parameterised for p=2 and ADER-O2/O3, plus a 200-step bulk stability sweep.

```
=== p=2 mixed-flux Adjacent reproducer ===

[Phase 1] Single-step Gate 2-analytic-ADER at p=1, p=2 with various ADER-O

----- p1_O2 (order=1, ader_order=2) -----
  ndof_total = 192  |central_flux_face_set_| = 32
  rel residual vs analytic = 3.976e-13       PASSED

----- p2_O2 (order=2, ader_order=2) -----
  ndof_total = 480  |central_flux_face_set_| = 32
  rel residual vs analytic = 1.518e-13       PASSED

----- p2_O3 (order=2, ader_order=3) -----
  ndof_total = 480  |central_flux_face_set_| = 32
  rel residual vs analytic = 1.518e-13       PASSED

[Phase 2] Multi-step stability (200 ADER steps, dt = 1e-5)
  [p1_O2_none]   max|Q|=4.677e+08  growth=1.218
  [p1_O2_adj]    max|Q|=4.854e+08  growth=1.264
  [p2_O2_none]   max|Q|=4.232e+08  growth=1.102
  [p2_O2_adj]    max|Q|=4.166e+08  growth=1.085
  [p2_O3_none]   max|Q|=4.232e+08  growth=1.102
  [p2_O3_adj]    max|Q|=4.166e+08  growth=1.085

[Phase 3] Stability ratios — Adjacent / None
  p1_O2 ratio = 1.038
  p2_O2 ratio = 0.984
  p2_O3 ratio = 0.984
```

**Verdict from the reproducer:**
- The single-step **algebraic identity** `Q_new_adj − Q_new_none == M⁻¹ Σ_f w · shape · (F_int(I) − F_central(I))` holds to ~1.5e‑13 relative at **all** tested combinations (p=1/p=2, O2/O3). The Mixed-Flux dispatch IMPLEMENTATION is mathematically correct at p=2.
- The **bulk-only** 200-step stability sweep is comparable across p=1/p=2 and Adjacent/None modes — no blow-up on the bulk wave operator alone.

These two pieces of evidence rule out a code bug in the bulk Mixed-Flux dispatch as the cause of the production blow-up.

---

## Implementation audit (no defects found in the bulk dispatch)

| Component | File:line | Status |
|-----------|-----------|--------|
| Central flux primitive `F = ½ T·(Aₓ⁺+Aₓ⁻)·Tinv·(Q_self+Q_nbr)` | `dynamic/godunov_flux.cpp:382-423` | Correct (verified by `test_godunov_central_flux`) |
| `BuildCentralFluxFaceSet_` — Adjacent walks fault-adjacent elements, inserts non-fault interior+shared faces | `dynamic/wave_operator.inl:1411-1565` | Correct (verified by `test_mixed_flux_face_set`) |
| Mult/RK1 dispatch site | `dynamic/wave_operator.inl:2642-2650` | Identity verified |
| ADER local-face dispatch | `dynamic/wave_operator.inl:4158-4172` | Identity verified at p=1, **p=2, O3** |
| ADER shared-face dispatch | `dynamic/wave_operator.inl:4672-4691` | Same kernel; structurally identical |
| MPI shared-face exchange (Allgatherv keys) | `dynamic/wave_operator.inl:1470-1565` | Correct |
| Driver → CFL normalisation `cfl_factor / (3·(2p+1))` | `drivers/tpv*_driver.cpp:1059` | Correct (p-order included) |

The audit above included reviewing every dispatch site, the rotation pipeline, the per-DOF Jacobian application in the CK recursion, and the symmetric face-quadrature integration at face order `2*order_`.

---

## Time integrator: SEAS-MFEM vs drdg3d (reference)

The user asked to compare. Findings:

| Aspect | SEAS-MFEM | drdg3d (`mod_rk.F90`) |
|--------|-----------|------------------------|
| Time integrator | **ADER-DG** with Cauchy-Kovalevskaya predictor: `I = Σ_{k=0}^{O−1} (dt^{k+1}/(k+1)!) D(k)` where `D(k+1) = −Aᵈ ∂_d D(k)` (`wave_operator.inl:1004-1083`) + corrector | **Explicit Runge-Kutta** (4-stage RK4 / 5-stage LSRK / SSPRK3 / SSPRK2; `mesh%nrk = {2,3,4,5}` selected by `timeIntegrationMethod`) |
| Predictor coupling | **Element-local** — predictor sees no neighbour state | Each RK stage is fully coupled (Q on every element gets the latest neighbour fluxes per stage) |
| Central-flux usage at fault-adjacent face | `flux_.Central(nor, I_self, I_nbr, F_h)` in corrector (uses time-integrated state) | `fstar` from `dT_x+dT_y+dT_z`-style global jump (`mod_wave.F90:587-602`) inside an RK stage on the **instantaneous** state |
| Higher-p extension | Untested at p≥2 with mixed flux (driver guards against it; see §Production guard) | Validated up to `Order = 3` (cubic) |

**Why this matters for the bug:** drdg3d's central flux is evaluated on the instantaneous Q at every RK stage. Each stage has its own face evaluation, providing many opportunities to dissipate high-frequency modes via the upwind portion of the operator. SEAS-MFEM's ADER predictor is element-local — the time-integrated state `I` carries `dt² L(Q)` and (at p=2) `dt³ L²(Q)` terms whose neighbour coupling is ONLY restored by the corrector's face flux. With central flux at adjacent faces, that single face-coupling pass at the corrector lacks the dissipation drdg3d's RK gets across multiple stages.

This is a structural difference, not a code bug, but it explains why drdg3d's mixed-flux works at higher p while SEAS-MFEM's exact same dispatch does not.

---

## Root cause (most likely)

The bug is **NOT in the bulk Mixed-Flux dispatch** — that part is mathematically correct (verified) and bulk-stable (verified). The catastrophic blow-up arises from the **interaction between Mixed-Flux Adjacent dispatch and the friction-coupled fault flux at p≥2**.

Concretely:

1. **High-frequency modes accumulate in fault-adjacent elements at p=2.** With central flux on the faces between fault-adjacent elements and their non-fault neighbours, the operator there is non-dissipative. At p=2 the basis supports more high-frequency modes (10 DOFs per tet vs 4 at p=1) which the central-flux corridor cannot damp.

2. **The fault face's per-side imposed-state flux** (`fault_flux_->EvaluateADER` at `wave_operator.inl:3627`) **consumes a trial state derived from the local I**:
   ```
   sigma_n_trial = η_p · (Q⁻[VX] − Q⁺[VX] + σ_n⁺/Zp⁺ + σ_n⁻/Zp⁻)   (eq 7a)
   τ_1_trial    = η_s · (Q⁻[VY] − Q⁺[VY] + τ_1⁺/Zs⁺ + τ_1⁻/Zs⁻)
   τ_2_trial    = η_s · (Q⁻[VZ] − Q⁺[VZ] + τ_2⁺/Zs⁺ + τ_2⁻/Zs⁻)
   ```
   `fault_face_flux.cpp:53-66`. When the trial state is amplified by accumulated high-frequency modes from the central corridor, the friction solver at fault QPs returns a runaway slip rate.

3. **Driver guard already flags this as untested**: `drivers/tpv*_driver.cpp:602-628` (R-1204) FATAL-aborts on `--mixed-flux ≠ none` + `--ader-order > 2` unless `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1`. The failing sbatch sets this override. Comment cites: *"Mixed-flux dispatch was validated at ADER-O2 only."*

4. **`ComputeMaxDt`'s 0.9× factor for Adjacent is too generous** (`wave_operator.inl:4983-5001`). Zhang 2023 §3.3 cites CFL ≈ 0.3 for mixed-flux runs; the current 0.9 factor only reduces dt by 10%. The factor was set based on bulk-only stability arguments and does not account for fault-coupling amplification at higher p.

---

## Hypotheses ranked

| # | Hypothesis | Status | Test that confirms / rules out |
|---|-----------|--------|--------------------------------|
| 1 | Bug in Central() flux primitive | **RULED OUT** | `test_godunov_central_flux` (Gates 1–3 pass) |
| 2 | Bug in central_flux_face_set_ construction (wrong faces marked central) | **RULED OUT** | `test_mixed_flux_face_set`, p=2 reproducer single-step identity |
| 3 | Bug in mode-dependent dispatch logic at face loop | **RULED OUT** at p=2 | New `test_mixed_flux_dispatch_adjacent_p2` Phase 1: identity holds to 1.5e‑13 |
| 4 | Bulk wave operator unstable at p=2 + Adjacent | **RULED OUT** | New reproducer Phase 2: 200 steps, ratio Adjacent/None ≤ 1.04 |
| 5 | Friction-coupled fault dynamics + central-flux corridor unstable at p≥2 | **LIKELY** | Production runs blow up; bulk-only is stable; consistent with §Root cause |
| 6 | CFL factor `0.9` for Adjacent is too loose at p≥2 with friction | **POSSIBLE CONTRIBUTOR** | Tightening factor to ~0.3 + retesting would confirm/refute |

---

## Recommended fixes (ranked, lowest risk first)

### Fix A (immediate gate hardening — zero numerical risk)

Tighten `wave_operator.inl:4983-5001` to scale the Adjacent factor with polynomial order when the fault is engaged. Concretely:

```cpp
case MixedFluxMode::Adjacent:
   // R-NEW: empirical safety factor for friction-coupled runs.  For
   // p≥2 the central-flux corridor near the fault accumulates
   // high-frequency modes that the per-side imposed-state Riemann
   // cannot fully damp; tighten CFL to compensate.  See debug doc
   // mixed_flux_p2_debug.md §Root cause.
   cfl_mixed_flux_factor = (order_ >= 2) ? 0.3 : 0.9;
   break;
```

### Fix B (re-enable the FATAL guard)

Remove the `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2` override path entirely (or change it to a hard FATAL with no override) until a fault-coupled stability gate is added. The override silently lets users run an unverified scheme.

### Fix C (add a fault-coupled stability gate)

Build a TPV205-style (LSW friction, simplest fault coupling) smoke test with:
- 4–48 tet mesh with embedded fault (BuildSmallFaultTetMesh fixture is suitable)
- LSW friction parameters from TPV205 (`config/tpv205_params.hpp`)
- 200 ADER steps with p∈{1,2}, mode∈{None, Adjacent}, ADER-O∈{2,3}
- Pass criterion: max-slip-rate growth ratio Adjacent/None ≤ 2× across all configurations

If this gate fails at p=2 + Adjacent with friction (which the production runs strongly suggest), then the fix path is to either tighten CFL (Fix A), add a partial-upwind blend at adjacent faces, or restrict mixed-flux to p=1 in production.

### Fix D (mathematical fix — α-blend)

Zhang 2023 also discusses a weighted central-upwind blend. Replace the binary Central/Interior switch at adjacent faces with `F = α · F_int + (1−α) · F_central` where α is a small upwind weight (e.g. α = 0.1 at p=2). This preserves most of the Mixed-Flux dispersion benefit while restoring O(α) dissipation. Requires a literature/scheme review and a stability calibration; do not adopt without an analysis pass.

---

## Files modified by this debug session

- **NEW** `tests/unit/test_mixed_flux_dispatch_adjacent_p2.cpp` — reproducer.
- **NEW** `debug_document/mixed_flux_p2_debug.md` — this document.

No production source files were modified. The reproducer is intentionally isolated from production code.

---

## Verification of conclusions

The reproducer test was built with `mfem-dev` conda env (`mpicxx` from `/Users/chunhuizhao/miniforge/envs/mfem-dev`) and run locally:

```
./seas_test_mixed_flux_dispatch_adjacent_p2
...
Results: 3 passed, 0 failed out of 3 tests
```

All three Gate 2-analytic-ADER identity checks PASS at p=1 (O2), p=2 (O2), and p=2 (O3). The bulk multi-step sweep produced no blow-up at any configuration.

The conclusion that the bulk Mixed-Flux dispatch is correct at p=2 is therefore directly evidence-backed, not inferred.
