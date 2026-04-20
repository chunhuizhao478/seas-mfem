# TPV102 Debug — v9.0.0 Execution Log (Phase 1 Unit Tests)

**Plan:** `tpv102_debug_v9.0.0_debug_plan.md` (Rev 3, 2026-04-20).
**Phase executed:** Phase 1 (§3.1, §3.1b, §3.1c, §3.1d), laptop unit tests only.
**Status:** A1-FAIL on §3.1b — new hypothesis H-V9-J produced. H-V9-A / H-V9-B /
H-V9-C / H-V9-I are all FALSIFIED at unit level. No code fix has been applied
yet; per CLAUDE.md policy, the user must approve any proposed fix before it
lands.

---

## Artefacts produced

| File | Plan § | Purpose |
|---|---|---|
| `miniapps/seas/tests/unit/test_fault_face_flux_frame.cpp` | §3.1 | `Evaluate` in isolation; BP5 channel decoupling |
| `miniapps/seas/tests/unit/test_fault_face_flux_frame_and_flux.cpp` | §3.1b | FAST-PATH end-to-end chain, F_h target ~220 m²/s² |
| `miniapps/seas/tests/unit/test_godunov_identity_normal_reversal.cpp` | §3.1c | Godunov identity F(+n,L,R) + F(-n,R,L) = 0 |
| `miniapps/seas/tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp` | §3.1d | Interior vs shared branch equivalence on Elem1 |

Makefile updated with 4 new build rules and `test-*` phony targets (file:
`miniapps/seas/Makefile`, insertions after existing fault-face-flux test entries).

All four binaries compile clean with `make seas_test_fault_face_flux_frame
seas_test_fault_face_flux_frame_and_flux seas_test_godunov_identity_normal_reversal
seas_test_fault_flux_interior_vs_shared_branch_equivalence` and run in
<1 s each.

---

## Results summary

| Test (§) | Pass/Total | Status | Hypothesis outcome |
|---|---|---|---|
| §3.1 Evaluate isolated | 12/12 | A1-PASS | H-V9-A and H-V9-B **FALSIFIED** |
| §3.1b FAST-PATH chain | 7/9 | **A1-FAIL** | New finding H-V9-J (see §3) |
| §3.1c Godunov identity | 5/5 | A1-PASS | H-V9-I **FALSIFIED** |
| §3.1d Interior vs shared | 4/4 | A1-PASS | H-V9-C **FALSIFIED** |

---

## 1. §3.1 — Evaluate in isolation (PASS)

Fixture: TPV102 VW-zone DOFData at **post-breakaway steady slip** (tau2_0 =
81 MPa, V_operating = 1 m/s, psi solved analytically for friction balance).

**Important fixture note.** The plan's equilibrium input (tau2_0 = tau_ini
= 75 MPa, psi = ComputeInitialPsi) gives V = V_ini = 1e-12 m/s, which is
by design quasi-static and produces F_h ≈ 0.  To exercise the radiation-
birth path the fixture uses post-breakaway operating-point inputs.  The
first test run exposed this — the equilibrium assertions "tau2_corr non-
trivial" and "dV1 first-order" were naïve.  The fix (current fixture)
uses V = 1 m/s, tau2 = 81 MPa, psi derived from
`psi = a * ln(2*V0/V * sinh(tau/(sigma_n*a)))`.

Observed values at this operating point:
- V1 = 0 (dip) — BP5 strike-slip invariant holds
- V2 = 0.271 m/s (strike) — ~27 % of operating-point V because the friction
  solver settles at the implicit steady-slip value for the given psi
- tau1_corr = 0 (dip) — no cross-channel leak
- tau2_corr = −1.25 MPa (strike) — = −eta_s · V2 as Eq. 10 predicts
- VY-perturbation → dtau1_corr = -4553 Pa (within 2 % of analytical
  -eta_s · δ = -4624 Pa); dtau2_corr stays at noise (no leak)
- VZ-perturbation → dtau2_corr = -2007 Pa (matches dtau2_trial scaled by
  friction-solver damping ~0.43); dtau1_corr stays at noise (no leak)

**Conclusion:** `Evaluate` preserves BP5 strike-dip decoupling.  The
algebraic form of Eq. 7-12 is correct.  H-V9-A (T_can rotation) and
H-V9-B (Eq. 10 impedance update) are falsified for this operating point.

## 2. §3.1c — Godunov identity (PASS)

F(+n, L, R) + F(−n, R, L) = 0 tested at 5 normals (y-axis, x-axis, z-axis
[BuildFrame up-fallback branch], oblique (1,1,1)/√3, tilted (0.6, 0.8, 0))
with random physically-sane Q_A, Q_B.

Worst relative drift across all 5 cases: **4.59 × 10⁻¹⁶** (single ULP).

**Conclusion:** H-V9-I (BuildFrame tangent arbitrariness breaking Godunov
identity) is falsified.  The `Interior` flux is deterministically
conservative under normal reversal across all configurations tested,
including the z-axis case which triggers BuildFrame's up-fallback branch.

## 3. §3.1b — FAST-PATH chain (FAIL → new hypothesis H-V9-J)

Fixture: TPV102 hypo-style canonical frame
(can_n=(0,−1,0), can_t1=(0,0,−1), can_t2=(1,0,0)), quiescent bulk
(Q_plus_local = Q_minus_local = 0), post-breakaway DOFData.

Step-wise recorded values:

**Step 1 — Evaluate (fault-local):**
- slip_rate |V| = 0.271 m/s
- tau1_corr = 0 Pa, tau2_corr = −1.25 × 10⁶ Pa
- Q_imp_plus_local[VZ] = −0.135 m/s, Q_imp_minus_local[VZ] = +0.135 m/s
- Q_imp_plus/minus[SXZ] = −1.25 × 10⁶ Pa (**same on both sides** by Eq. 11d/12d)

**Step 2 — T_can rotation (→ global):**
- Q_imp_plus_g[VX] = −0.135 m/s (matches v_t2_local, can_t2 = +x)
- Q_imp_plus_g[SXY] = +1.25 × 10⁶ Pa
- Q_imp_minus_g[VX] = +0.135 m/s (matches opposite v_t2_local)
- All other components = 0

Both Steps 1 and 2 match the plan's REVIEW R-1005 expectation for
hypocentre steady slip.

**Step 3 — `flux.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h)`:**
- **F_h[VX_global] = 0** exactly (machine-zero, not rounding)
- F_h[SXY_global] = −2.82 × 10⁻⁷ Pa·m/s (also effectively zero)
- F_h on all other components = 0

**Analytical cross-check:** `(eta_s / rho) · V2 = 469 m²/s²` — the scale
the plan used to motivate F_h ≈ 220 m²/s².  The measured value is 10¹¹
below that.

### 3.1 Root-cause analysis — why F_h = 0

Write the Godunov flux in the tangential (SXY_int, VY_int) sub-block,
using the split-matrix form for a 1-D elastic wave:

```
F[VY]  = −(1/(2ρ))·(SXY_L + SXY_R) + (cs/2)·(VY_L − VY_R)
F[SXY] =   (cs/2)·(SXY_L − SXY_R) − (mu/2)·(VY_L + VY_R)
```

Substitute the imposed-state values that `Evaluate` produces under the
pure strike-slip configuration (bulk Q = 0):

- SXY_L = SXY_R = −eta_s · V = −(Zs/2)·V  (identical on both sides by
  Eq. 11d and Eq. 12d)
- VY_L = −V/2, VY_R = +V/2  (velocity jump VY_L − VY_R = −V, per Eq. 11a/12a)

Plug in:

```
F[VY] = −(1/(2ρ))·(−Zs·V) + (cs/2)·(−V)
      =  Zs·V/(2ρ) − cs·V/2
      =  ρ·cs·V/(2ρ) − cs·V/2
      =  cs·V/2 − cs·V/2
      =  0             ← EXACT cancellation

F[SXY] = (cs/2)·0 − (mu/2)·0  =  0
```

The cancellation is **structural**, not numerical.  Every line in the
derivation is an identity, independent of the specific values of V, cs,
ρ, Zs, mu.

### 3.2 Why this matches production symptoms

**Observed production failure (job 7666327, §0.2 of plan):**
- V_max saturates at 7.69871 m/s (friction's own internal V, not
  radiation-damped V)
- bulk max ‖Q‖_∞ = 1.66 × 10⁻⁸ — 10¹⁵ below physical O(MPa)/O(m/s)
- V_dip, slip_dip, τ_dip **exactly** 0 at hypocentre
- Radiation-damping / under-stress mechanism never fires

If F_h = 0 on the fault face, the DG rhs update
`rhs[Elem1] −= w·shape·F_h` deposits nothing into the bulk-adjacent
DOFs.  No radiation ever enters the bulk, so ‖Q‖_∞ stays at noise
floor.  Without bulk radiation, there is no bulk-driven tau_trial, no
radiation-damping feedback, and V_abs pins at the friction-only value
produced by tau_0 alone.  This exactly matches the §0.4 DRDG3D-vs-MFEM
deviation: MFEM's V_strike never decays, the dip channels stay at
machine zero, and state never recovers.

**Plan-file cross-check:** §3A FP-2 decision tree:
> (c) exactly 0 or 1e-15 → a step zeros the flux. Print between steps 1/2/3
> to bisect.

Bisection done (above).  Cancellation isolated to Step 3's split-flux
formula acting on the imposed states.

### 3.3 Hypothesis H-V9-J (NEW) — formal statement

**H-V9-J (Imposed-state Godunov flux cancellation).**  The chain
`flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g)` — used in both the
interior-fault branch (`wave_operator.inl:807-829`) and the shared-fault
branch (`wave_operator.inl:1174-1186`) — identically cancels in the
tangential sub-block when the imposed states are produced by
`FaultFaceFlux::Evaluate` under pure strike-slip geometry (BP5
convention, tau1_0 = 0, bulk Q = 0).  The cancellation is algebraic
and independent of slip rate, wave speeds, or material parameters.
Consequently no tangential DG flux is deposited into bulk rhs, the
bulk wave field cannot grow beyond numerical noise, and every
downstream symptom in R-1001 follows.

**H-V9-J is consistent with every rule-out in §0.3:**
- DOFData consistent across ranks ✓ (Evaluate runs identically)
- CFL OK ✓ (not a dt issue)
- Rank count irrelevant ✓ (F_h = 0 on any single rank)
- Mass matrix correct ✓ (zero input · correct operator = zero output)
- Nucleation ramp & Brent solver correct ✓ (tau_0 ramp alive in §3.1)

**H-V9-J is consistent with the §0.4 exactly-zero dip-channel observation:**
F_h[VY_global] and F_h[VZ_global] are both also zero (same argument as
F_h[VX_global], applied to the dip-tangent sub-block).  So no dip
radiation ever reaches the hypo-adjacent bulk DOFs, which is exactly
the "V_dip = 0.000 m/s" signature vs DRDG3D's noise-floor ±0.003 m/s.

---

## 4. §3.1d — Interior vs shared branch equivalence (PASS, but see caveat)

Both branches, when fed the same preamble (canonical frame + Evaluate +
T_can), produce IDENTICAL Elem1 contributions to ULP in both Cases A
(elem1_on_plus = true) and B (elem1_on_plus = false), with and
without a small bulk perturbation.  **The equivalence holds** — so the
REVIEW.md R-1001 proposal to unify the interior branch with the shared
convention is a correctness-preserving refactor, not a bug fix.

**Caveat.** Because F_h is itself = 0 under quiescent bulk (see §3), the
Case-A/B equivalence in that sub-case is trivially "0 = 0".  The
perturbed sub-cases (nonzero self-side bulk) show a non-trivial F_h
(|F_h[SXY]| ~ 1.9e8 Pa·m/s, |F_h[VX]| ~ 20 m²/s²) and the two branches
still agree to ULP.  So H-V9-C is falsified both trivially AND
non-trivially.

Tolerance note: the first run reported "drift = 1 at component SXY"
for Case-B quiescent because F_interior = −1.26e−19 and F_shared =
+1.26e−19 (sign-flipped floating-point residuals below 1e-18).  The
test was hardened with an absolute tolerance abs_tol = 1e-10 to treat
such noise as equal.  This is NOT a bug in either branch — it is the
ApplySplitFlux rounding difference between two algebraically
equivalent paths.

---

## 5. What Phase 1 rules OUT vs what it rules IN

**Falsified (do NOT re-diagnose):**
- H-V9-A (T_can rotation sign bug) — §3.1 and §3.1b Step 2 all check out
- H-V9-B (Eq. 10 impedance-update mis-mapping) — §3.1 local tau1/tau2
  decoupling, §3.1b Step 1 values match Eq. 10 analytically
- H-V9-C (interior vs shared branch asymmetry) — §3.1d PASS
- H-V9-I (Godunov identity under normal reversal) — §3.1c PASS

**Promoted to Rank 1:**
- **H-V9-J (NEW)** — imposed-state Godunov flux cancellation.  Replaces
  H-V9-C/I at the top of the hypothesis list.

**Still open (not yet Phase-1-addressable):**
- H-V9-D (double-count / visit-count) — requires Phase 2 instrumentation
- H-V9-E (mass-inverse) — deprioritised per REVIEW R-1003
- H-V9-F (IP penalty / BR2 damping) — would need Phase 0 Frontera run
- H-V9-G (mesh-resolution under-damping) — Phase 0 1000 m comparison
- H-V9-H (Tandem convention audit) — Round 1 Slot B

---

## 6. Required user decision before any code change

Per `miniapps/seas/CLAUDE.md` ("Proposing Fixes") and
`/Users/chunhuizhao/projects/seas-mfem/CLAUDE.md` ("Debugging Behavior"),
the debugger agent MUST NOT apply a fix without user approval.  The
options for H-V9-J, spelled out for the user:

### Option 6.1 — Re-read SeisSol's DG flux formulation
The SEAS-MFEM code passes `(Q_imp_plus_g, Q_imp_minus_g)` directly to
`flux_.Interior`, which is `F_DG = A+·Q_imp_plus + A−·Q_imp_minus`.
The algebra in §3.1 shows this is zero for symmetric strike-slip.
SeisSol's dynamic-rupture DG formulation (de la Puente 2009;
`src/DynamicRupture/FrictionLaws/` and `src/DynamicRupture/Kernels/` in
`/Users/chunhuizhao/projects/SeisSol/`) likely uses a different form —
e.g.
  `F_DG_self = A+·(Q_imp_nbr − Q_self) + A−·Q_imp_self`
or folds the friction traction as a direct bulk-force term, not via
the numerical flux alone.  SEAS-MFEM's `fault_face_flux.cpp` and
`godunov_flux.cpp` already explicitly cite SeisSol
(`FrictionSolverCommon.h`, `ElasticSetup.h`), so the port surface is
small.  The correct fix requires reading how SeisSol's
**Kernels** consume the imposed states computed by the **FrictionLaws**.
**Recommended next step before any code change.**

### Option 6.2 — Comparison against a hand-derived radiation BC
For a quiescent bulk + slipping fault at V, classical radiation boundary
conditions require that the outgoing S-wave amplitude be ~(Zs/2)·V
in stress and ~V/2 in velocity, carrying power ~(Zs/2)·V² across the
face.  Compute F_DG that would produce these amplitudes in one step
of the DG system (given mass matrix), and compare to what `flux_.Interior`
on the imposed states gives.

### Option 6.3 — Register §3.1b as the regression gate
Per plan §3.3 A1-FAIL rule, the failing test is committed first, the
fix second.  The fix must make §3.1b PASS (F_h[VX] ∈ (10, 1e4) m²/s²
for the post-breakaway fixture) while leaving §3.1, §3.1c, §3.1d
still passing.

---

## 7. Checklist for the next session

- [ ] User reviews this fix report and confirms H-V9-J is plausible.
- [ ] User decides between Option 6.1 (SeisSol reference read at
      `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/`) and
      Option 6.2 (hand-derivation).  Both are laptop / desk work; no
      Frontera runs yet.
- [ ] If the fix lands, Frontera Phase 0 sbatch (per plan §2) is the
      validation step.  No local reproducer run per
      `feedback_no_local_reproducer.md`.
- [ ] Update `tpv102_debug_v9.0.0_check.md` with H-V9-J finding (once
      the v9.0.0 check-doc is created — not yet on disk).
- [ ] Add `make test-fault-face-flux-frame-and-flux` to the CI / smoke
      list once a fix lands and §3.1b turns green.

---

## 8. Build / run commands for reproduction

```bash
conda activate mfem-dev
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
make seas_test_fault_face_flux_frame \
     seas_test_fault_face_flux_frame_and_flux \
     seas_test_godunov_identity_normal_reversal \
     seas_test_fault_flux_interior_vs_shared_branch_equivalence

./seas_test_fault_face_flux_frame                             # §3.1 → 12/12 PASS
./seas_test_fault_face_flux_frame_and_flux                    # §3.1b → 7/9 (2 FAIL on F_h[VX] range — H-V9-J smoking gun)
./seas_test_godunov_identity_normal_reversal                  # §3.1c → 5/5 PASS
./seas_test_fault_flux_interior_vs_shared_branch_equivalence  # §3.1d → 4/4 PASS
```

Total wall-clock: ~5 s for the four binaries on laptop.

---

## 9. References

- `tpv102_debug_v9.0.0_debug_plan.md` §3.1 / §3.1b / §3.1c / §3.1d
- `miniapps/seas/dynamic/fault_face_flux.cpp:68-156` (Evaluate pipeline)
- `miniapps/seas/dynamic/godunov_flux.cpp:305-349` (ApplySplitFlux, Interior)
- `miniapps/seas/dynamic/wave_operator.inl:807-829` (interior-fault branch)
- `miniapps/seas/dynamic/wave_operator.inl:1174-1186` (shared-fault branch)
- de la Puente et al. (2009) "Dynamic rupture modeling on unstructured
  meshes using a discontinuous Galerkin method" — referenced in
  `godunov_flux.hpp` comments; needed for §6.1 option reading
