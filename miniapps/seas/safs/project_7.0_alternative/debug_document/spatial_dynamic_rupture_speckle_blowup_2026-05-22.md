# Spatial dynamic rupture — on-patch "speckle" blow-up (job 7744292, 2026-05-22)

Second SAFS dynamic-rupture blow-up, **distinct** from the 2026-05-20
north-tip sliver event (`spatial_dynamic_rupture_postevent_debug_2026-05-20.md`).
Run: `spatial_dyn_smoke` job 7744292, 400 ranks, config
`config/spatial_friction_slip_weakening_safs_projected_stress.toml`
(psi=37° / az351 H&Z stress, 500 m `z0embed` sliver-free mesh, μ_s=0.65 /
μ_d=0.30 / D_c=1.0, gradual_overstress δτ=6.5 MPa).

---

## Symptom

`slip_rate_strike` on the fault shows a **smooth Gaussian nucleation patch
(correct)** surrounded by **salt-and-pepper** — isolated, grid-scale DOFs
with huge slip rate, NOT a coherent propagating front (ParaView, t≈0.84 s).
V_max trace:

| step | t (s) | V_max (m/s) | regime |
|---|---|---|---|
| 1000–2300 | 0.35–0.81 | ~3 | healthy nucleation patch slipping |
| 2400 | 0.84 | 975 | blow-up onset |
| 2400–12900 | 0.84–4.51 | 975 → 5.7e145 | exponential growth ~3.2%/step |
| 13000–13200 | 4.55–4.62 | 5.4e136 → 1.5e89 → **3.17** | **recovery** |

Key features distinguishing this from the 2026-05-20 event:
1. **Early** (t≈0.8 s, during the nucleation regime) — not post-event (t≈3 s).
2. **On / next to the nucleation patch** — not at the far north tip.
3. **Recovers** to the physical ~3 m/s — the 05-20 event went to NaN.
4. Mesh is the new sliver-free `500m z0embed`, so it is **not** the slivers.

The recovery (finite the whole way, 1e145 → 3, never NaN) rules out a global
CFL / scheme instability (those never recover). It is a **localized,
self-limiting** instability: once the speckled DOFs slip past D_c = 1.0 m,
μ pins to μ_d, the slip-weakening gradient that feeds the growth is
exhausted, the corrupted stress radiates away, and the genuine rupture front
(~3 m/s, elsewhere) becomes the global max again.

---

## Evidence — the dynamic-rupture CODE is correct (TPV205 comparison)

Traced both drivers. The fault machinery is the **same code path** for the
proven planar TPV205 driver and the SAFS spatial driver:

`Tpv205SubStepIterator::AdvanceWithSubStepStates` → `StepOneQP_` →
`FaultFaceFlux::ComputeTrialTraction` / `SolveLSW_TPV205` /
`BuildImposedState`, with shared-fault QPs through
`FaultFaceFlux::EvaluateADER_LSW`, corrector via `WaveOperator::AdvanceADER`.

- `SolveLSW_TPV205` (tpv205_friction.hpp:93) is the standard radiation-damped
  balance `V = max(0, (|τ|−μσ_n)/η_s)`, V ∥ τ — stable, no sign error.
- Slip accumulation `d.slip{1,2} += V*dt_sub` (tpv205_substep_iterator.cpp:122)
  is correct and single-owner (the iterator; `EvaluateADER_LSW` deliberately
  does NOT re-accumulate — fault_face_flux.cpp:798).
- CFL is reduced by the identical `cfl/(3·(2N+1))` factor as TPV205
  (spatial_dyn_driver.cpp:1218; the job-7743554 lesson). dt is NOT the issue.
- Interior flux is plain homogeneous Godunov (`--no-sidecar-material`,
  λ=μ=32 GPa); the Phase-R bi-material path is inactive.
- Initial stress is smooth in ParaView before rupture → the per-DOF stress
  projection (`FaultGeometry::ComputeSAFSParams`) and the gradual_overstress
  resolver are NOT injecting salt-and-pepper into the IC.

**No bug found in the dynamic-rupture code.** TPV205 is stable on the same
code because of three protections SAFS removed:

| Protection | TPV205 (stable) | SAFS run 7744292 (blows up) |
|---|---|---|
| Cohesive-zone resolution | well-resolved | `L_nuc/h_min` min = **4.79** (criterion ≥ 10) |
| Strength barrier outside rupture area | `mu_s_barrier`=1e4 → `V≡0` | none on the active fault (only z∈[−20,−15] km) |
| Background criticality | uniform, subcritical | `\|τ_pre\|/(μ_s·σ_n_eff)` outside = **0.986** |

---

## Root cause

A **two-factor** numerical instability (same class as the 05-20 doc, but
moved early and onto the patch by the new psi37 calibration):

1. **SEED** — the cohesive zone is under-resolved (`L_nuc/h_min` min 4.79 ≪ 10;
   ~2.5 elems in the patch). An under-resolved slip-weakening front emits
   grid-scale (salt-and-pepper) slip-rate oscillations of a few %.
2. **AMPLIFIER** — the background is 1.4% below static yield (outside ratio
   0.986) with no barrier to clamp it. The few-% grid noise tips near-critical
   background DOFs over yield; slip-weakening positive feedback then grows
   each one exponentially (~3.2%/step) until D_c saturation halts it.

The `[derived]` well-posedness gate (`spatial_print_derived.cpp`) gates only
on `outside_max ≥ 1.0` (0.986 squeaked by) and `L_nuc/h_min` is **printed,
never gated** — so the run launched straight into the known blow-up regime.

### Spatial evolution (ParaView, slip_rate_strike) — it is a TRAVELING band
Frames at t = 0.94 / 1.64 / 2.09 / 3.23 / 4.57 s show the checkerboard is NOT
static — it RIDES the rupture front:
- t=0.94: ignites **down-dip of the patch** (smooth Gaussian intact).
- t=1.64: semicircular checkerboard on the hypocenter + a **detached** blob
  ahead along-strike (→ the near-critical background also ignites REMOTELY,
  not only at the front).
- t=2.09–3.23: a checkerboard band with a sharp front sweeps north along strike.
- t=4.57: the southern region has **HEALED to dark blue**; the band is mid-fault,
  front at the **geometric bend/stepover**.

Mechanism: the band ignites where the fault is actively weakening, grows on the
near-critical background, **heals behind** (DOFs slip past D_c → μ pins to μ_d →
the weakening gradient feeding it is gone → it radiates away), and the front
**arrests at the geometric bend** at t≈4.5 s. That arrest IS the V_max recovery
(steps 13000→13200): the front stops at the kink, the last band heals, the
global max drops back to the genuine ~3 m/s. So this is a *rupture-front* grid
instability, not a static IC defect.

### Why the fix works — the seismic S-ratio
S = (μ_s·σ_n_eff − |τ_pre|)/(|τ_pre| − μ_d·σ_n_eff) at the patch:
- OLD: (32.0 − 29.2)/(29.2 − 14.8) = 2.8/14.4 = **0.20** (violently supercritical
  → fast rupture, razor-thin dynamic cohesive zone → grid checkerboard).
- NEW: (41.9 − 29.2)/14.4 = 12.7/14.4 = **0.88**, stress drop UNCHANGED (μ_d=0.30).

Raising μ_s lifts S 0.20 → 0.88, still **below the Andrews threshold ~1.77** so
spontaneous rupture is preserved, but far more stable. It helps two ways:
(1) no remote background ignition (the 0.75-subcritical bet → kills the detached
blobs); (2) slower rupture → the *dynamic* cohesive zone is Lorentz-widened at
lower speed, partially offsetting the static L_nuc shrinkage from the wider
(μ_s−μ_d). Whether the FRONT checkerboard is fully killed or merely bounded is
the open question the re-run answers.

---

## Fix (applied) — lower background criticality, keep spontaneous rupture

Per run directive (chose: lower criticality; keep μ_d=0.30 strong stress drop;
ensure spontaneous rupture). Edited
`config/spatial_friction_slip_weakening_safs_projected_stress.toml`:

| param | old | new | why |
|---|---|---|---|
| `mu_s_default` | 0.65 | **0.85** | outside ratio 0.641/0.85 ≈ **0.75** (25% below yield) → removes the amplifier |
| `mu_d_default` | 0.30 | 0.30 | keep strong +14.4 MPa dynamic stress drop (propagation) |
| `delta_tau_strike_pa` | 6.5e6 | **20.0e6** | reach the higher static yield AND build a supercritical core > L_nuc |
| `d_c_default` | 1.0 | 1.0 | unchanged |

**Spontaneous-rupture verification** (patch σ_n_eff≈49.3 MPa, |τ_pre|≈29.2 MPa,
μ_bulk=32 GPa, Gaussian `F(r)=exp(−(r/radius)²)`, radius=3000 m):

- L_nuc = 32e9·1.0/((0.85−0.30)·49.3e6) = **1.18 km** (was 1.86 km; smaller → easier to exceed).
- yield gap at center = μ_s·σ_n_eff − |τ_pre| = 41.9 − 29.2 = 12.7 MPa.
- supercritical core radius r_c = 3000·√(−ln(12.7/20)) = **2.02 km**.
- margin **r_c/L_nuc ≈ 1.71** (> 1 required) → ruptures spontaneously, does not arrest at T_nuc.
- center overshoot = 29.2 + 20 − 41.9 = **+7.3 MPa** ≥ 0 → `[derived]` TRIGGER gate passes.
- in-patch dynamic ratio = 29.2/(0.30·49.3) = **1.97 > 1** → STRESS-DROP gate passes.
- δτ=20 MPa is below the 25 MPa "gross overdrive" ceiling.

All four `[derived]` gates green by hand; the gate aborts at startup
(`abort_on_failure=true`) if any prediction is off, so a bad recalibration
fails fast rather than wasting a multi-hour run.

---

## Next-run checklist (`--print-derived`)

Expect (within mesh/projection rounding):
- `max OUTSIDE ratio` ≈ **0.75** (no longer the 0.986 WARNING).
- `L_nuc / h_min` min ≈ **3.0** (max ≈ 9.7) — *worse* than before (see risk).
- `nucleation overshoot` ≈ **+7 MPa (sufficient)**.
- `in-patch min |τ_pre|/(μ_d·σ_n_eff)` ≈ **1.97 (> 1)**.
- `in-patch max |τ_pre|/(μ_s·σ_n_eff)` ≈ **0.70 (< 1, locked)**.

## Known remaining risk + next lever
This fix removes the **amplifier** (criticality), not the **seed**
(resolution) — raising μ_s shrank L_nuc, so `L_nuc/h_min` is now ~3.0
(worse). The bet is that a 25%-subcritical background can no longer grow the
grid noise. **If residual front speckle reappears** (watch `slip_rate_strike`
near the front in ParaView, and V_max for any super-linear growth), the next
lever is the **250 m z0embed mesh** (`[mesh].path`, ~2× finer →
`L_nuc/h_min` ~6–19), which attacks the seed directly. Pairing both is the
robust combination.

## Status
- CONFIRMED (code read): fault machinery identical to proven TPV205; LSW
  solve, slip accumulation, CFL all correct.
- APPLIED (config): μ_s 0.65→0.85, δτ_strike 6.5→20 MPa; spontaneous-rupture
  margin r_c/L_nuc ≈ 1.71 verified by hand.

---

## UPDATE — Round 2 (job 7745103): μ_s recalibration did NOT fix it

The μ_s=0.85 run (outside ratio confirmed 0.754, all `[derived]` gates PASS)
**still blows up**: V_max healthy ~4–6 m/s for ~1300 steps after nucleation,
then exponential growth from step ~2700 (t≈0.95 s): 58 → 1437 → 33066 m/s …
— essentially the same physical time as the μ_s=0.65 run.

**This ELIMINATES the criticality-amplification hypothesis.** A fault loaded
to only 0.754 of static yield cannot be tipped over by small radiated noise,
yet it blew up anyway. The Round-1 "near-critical background amplifier" story
was WRONG as the dominant cause. The growth is intrinsic to the actively-
rupturing region, independent of background stress level.

### New prime suspect: mixed-flux `adjacent` (non-dissipative central corridor)
- `debug_document/mixed_flux_p2_debug.md` ("p=2 Mixed-Flux Adjacent —
  Catastrophic Blow-up") documents THIS EXACT failure for TPV104/TPV205 +
  adjacent: slip rate → 1e10–1e11 m/s, "isolated high values across the
  propagated area, not just at peaks — consistent with an instability that
  grows wherever a wave passes a fault-adjacent (central-flux) face."  That
  is the SAFS speckle verbatim.
- Mechanism (that doc, Hyp 5 LIKELY): the central flux on fault-adjacent
  faces is NON-DISSIPATIVE (`godunov_flux.cpp:382` = ½·A·(Q_self+Q_nbr), no
  |A| term), and the element-local ADER corrector's single face-flux pass
  cannot damp the HF modes the way drdg3d's multi-stage RK does.
- `ComputeMaxDt` still uses `cfl_mixed_flux_factor = 0.9` for Adjacent
  (`wave_operator.inl:5186`) — only a 10% dt cut; the doc's recommended
  ~0.3 (Zhang 2023 §3.3) was never applied.
- SAFS extends the doc's p=2 finding to p=1: the under-resolved curvilinear
  front (`L_nuc/h_min` = 3.05) supplies the grid noise that the
  zero-dissipation central corridor then amplifies, even at p=1.
- TPV205 ran stably with adjacent because it is well-resolved (little front
  noise) — the flux is not inherently fatal; the SAFS noise source is.

### Discriminator (in flight): the no-mixed-flux run
`spatial_dyn_smoke_nomixedflux_*.sbatch` with `--mixed-flux none` (upwind
everywhere → fully dissipative).  The mixed_flux_p2 doc predicts `none` is
stable.
- `none` bounded → CONFIRMS mixed-flux adjacent.  Fix: `mixed_flux = "none"`
  for SAFS (simplest), or keep adjacent + tighten its CFL factor to ~0.3
  (doc Fix A; preserves the SSO benefit at higher cost).
- `none` also blows up → eliminates the flux → under-resolution/geometry;
  fix = 250 m z0embed mesh.

### Note on the Round-1 fix
The μ_s=0.85 / δτ=20 MPa recalibration is NOT reverted: a 0.75 subcritical
background is still desirable (no spontaneous background failure), and the
spontaneous-rupture margin is intact.  But it is NOT sufficient on its own.

---

## UPDATE — Round 3 (job 7745138): mixed-flux is ELIMINATED too

`--mixed-flux none` (log confirms `central-flux faces = 0`, full upwind
dissipation everywhere) **also blows up**: V_max healthy ~4–6 m/s through
step ~2600, then exponential from step 2700 (t≈1.05 s): 22 → 303 → 4035 m/s.
Same onset, same character.

So the Round-2 prediction ("none will be stable") was WRONG, and **mixed flux
is eliminated** as the cause. One real signal survives: `none` grew **slower**
(~2.6 %/step) than `adjacent` (~3.2 %/step) — the central flux was a genuine
*amplifier* (consistent with `mixed_flux_p2_debug.md`) but **not the source**.
The instability survives full upwind dissipation.

### Two hypotheses now falsified by runs
- μ_s / criticality (Round 2): 0.754 subcritical still blew up.
- mixed flux (Round 3): `none` (upwind) still blew up.

The blow-up is **invariant** to both → the source is a common factor of all
three runs: the **500 m mesh (L_nuc/h_min = 3.05, under-resolved)**, μ_d/d_c,
ADER-O2, or the curvilinear surface-rupturing geometry.  The screenshots show
a checkerboard **riding the rupture front** (igniting down-dip of the patch,
not at the free surface), and it onsets ~0.5 s after nucleation — consistent
with the *accelerating* front's cohesive zone Lorentz-contracting below the
grid as the rupture speeds up.  **Leading surviving hypothesis: cohesive-zone
under-resolution.**

Ruled out by code read this round:
- Nucleation is **fault-only** — `ApplyGradualOverstressIncrement`
  (`spatial_nucleation.cpp:174`) writes only `dof_data[i].tau{1,2}_nuc`
  (fault QP DOFData), consumed as fault traction; never a bulk-Q body force.
- Free surface is `natural_attrs={102}` and the fault reaches z=0, but the
  speckle ignites at depth (down-dip), not at the surface tip.

### Round-3 test (applied to config): D_c 1.0 → 2.0 m
Single-variable resolution test on the SAME 500 m mesh: D_c=2.0 raises the
patch `L_nuc/h_min` from 5.99 to **~12** (clears the ≥10 requirement;
global ~6.1–19.4).  D_c doubles L_nuc (1.18→2.36 km), so the nucleation
radius is enlarged 3000→4000 m to keep the supercritical core above L_nuc
(`r_c≈2.70 km`, margin 1.14) — otherwise the rupture would go subcritical and
arrest, confounding the test.  δτ stays 20 MPa.
- **bounded** → confirms under-resolution is the source; production path is
  the 250 m z0embed mesh at a *physical* D_c (2.0 m is inflated).
- **still blows up** → eliminates static cohesive-zone resolution too; the
  source is then the geometry (curvilinear/surface-rupturing front) or a
  resolution-independent code issue at the ADER LSW front — escalate to the
  250 m mesh and a per-element onset-locator.

### Expected `[derived]` for the D_c=2.0 / radius=4000 run
- `L_nuc/h_min` min/max ≈ **6.1 / 19.4** (patch ≈ 12; the printed *min* 6.1 is
  the far high-σ_n region, not the patch).
- `outside ratio` ≈ 0.754, `overshoot` ≈ +6.3 MPa, `in-patch dynamic` ≈ 1.96,
  `in-patch static` ≈ 0.70 — all unchanged (D_c/radius don't affect them).

---

## ROOT CAUSE CONFIRMED + FIXED (2026-05-22, local tilted-fault test)

The "resolution" hypotheses above were all **wrong**: the blow-up is invariant
to D_c (Dc=10 over-resolved still blew up), cfl, μ_s, mixed-flux, and
dip-prestress because the source is an **MPI/geometry orientation bug**, not a
front-resolution problem.

### The discriminator test
`tests/unit/test_rupture_tilted_fault_serial_vs_parallel.cpp` runs the SAME
2-tet TPV102 fault problem **serial vs np=2**, sweeping the fault-trace tilt
`θ` (background stress rotated with the fault so the physical loading is
θ-invariant), and asserts serial==parallel to 1e-10.  Nucleation is a uniform
fixed scalar in the local strike channel — the only way serial/parallel can
disagree is an orientation/side flip.

**Result (before fix):** serial==parallel to ~1e-15 for θ ≤ 60°, then a clean
**sign flip** (worst_rel = 2.0, component SYZ; max|V| 0.1164 → 0.1202) for
θ ≥ 75°.

### The mechanism (diagnostic-confirmed to the digit)
NOT `FaultBasis::ComputeOrientedFrame`'s `sign_flipped` (changing that rule is
a no-op — the cross-rank flip is symmetric under any line-canonicalization).
The bug is the **+/- side determination** in `wave_operator.inl`
(`interior_fault_elem1_on_plus_`, lines 459-465; `shared_fault_elem1_on_plus_`,
lines 519-525):

```
elem1_on_plus = ( (elem1_c - face_c) · ref_normal  <  0 )      // OLD (buggy)
```

It projects the centroid offset onto the **hardcoded `ref_normal=(0,-1,0)`**.
For an asymmetric tet the offset `(elem1_c - face_c)` has a tangential part; as
the fault tilts, `ref_normal`'s tangential component makes that term overtake
the true normal-offset term and the comparison flips.  For the test fixture the
margin is `kL·(sinθ/12 − cosθ/4)`, zero-crossing at `tanθ = 3 → θ = 71.57°`
(measured: θ=60 margin −52.8 → true; θ=75 margin +15.8 → FALSE).  Past the flip
**both** ranks report `elem1_on_plus = false`, breaking the "exactly one side is
+" invariant on the shared face → wrong ± role → the fixed-sign nucleation is
applied with the wrong sign on one side → blow-up.  Planar TPV/BP5 faults have
`n_phys = ±ref_normal` (`|dot| = 1`) so the margin is O(element-size) and never
flips — which is why TPV102/104/205 are stable on the same code.

### The fix
Project the centroid offset onto the **canonical face normal**, not
`ref_normal`, in BOTH side-determination blocks:

```
CalcOrtho(ftr->Face->Jacobian(), fn);                              // raw normal
if (FaultBasis::NormalNeedsFlipToCanonical(fn, ref_normal, 3)) { fn.Neg(); }
elem1_on_plus = ( (elem1_c - face_c) · fn  <  0 );                 // NEW
```

`NormalNeedsFlipToCanonical` (new shared static in `fault_basis.hpp`)
canonicalizes the physical normal **line** ±n to a rank/QP-independent
direction: `sign(n·ref_normal)`, with a largest-magnitude-component fallback in
the `|n·ref_normal| → 0` degenerate band (θ→90°).  On a shared face CalcOrtho
gives opposite raw normals on the two ranks; both map to the SAME canonical
direction → exactly one rank gets `elem1_on_plus=true` at every angle.
`ComputeOrientedFrame`'s `sign_flipped` now also routes through this helper
(single source of the rule).  For planar faults `canonical normal == ref_normal`
(up to a positive |J_F| scale) so the stored bool is unchanged — **byte-exact**.

### Verification (local, `conda activate mfem-dev`)
| test | before | after |
|------|--------|-------|
| **tilted-fault serial-vs-parallel** (the fix target) | FAIL worst_rel **2.0** @θ=75° | **PASS 7.4e-14** at all θ (0..90°) |
| sign-flipped-elem1-on-plus-truth-table | — | PASS 40/0 |
| rupture-multistep-serial-vs-parallel (planar) | — | PASS 7.2e-16 |
| fault-basis / dip-strike-symmetry / qp-orthonormality | — | PASS 219/0, 17/0, 1/0 (real TPV102 mesh) |
| shared-fault-role-consistency / dof-data-consistency | — | PASS 5/0, 4/0 |
| interior-fault-flux-path (+asym +rupture-regime) | — | PASS 6/0, 5/0, 3/0 |
| shared-fault-flux-path-analytic / godunov-interior-equal-sides | — | PASS, 135/0 |
| tpv102 locked-fault / absorbing-eq / pepper / ader-smoke | — | PASS 4/0, 1/0, 6/0, 4/0 |

**Confirmed pre-existing (byte-exact baseline match, NOT caused by this change):**
- `adjacent-triangle-fault-uniformity` (serial+parallel) — known **pepper-bug
  reproducer** (Frontera 7672897); stashed-baseline gives the IDENTICAL failure
  values (tau1_corr **2.236**, tau2_corr 1.744e-6, σ_n 2.049e-6) → planar fault,
  my change is a verified no-op here.  `first-step-audit` (Gate 16b) is the same
  planar fixture/family.
- `r101-shared-fault` — aborts on a missing `SetAbsorbingBackground` precondition
  (test-setup gap); side-role coverage is provided by the two passing
  shared-fault consistency tests above.
- `serial-parallel-consistency` — macOS **MUMPS Bus error** in the antiplane
  (BP1/BP2) `AntiplaneDomainOperator::Solve` path; not the 3D SAFS path.

### Remaining
- Re-run **Dc10 on Frontera** post-fix to confirm end-to-end (the step-0
  cross-rank check passed there because that partition's shared faces were at
  well-conditioned angles; the real fix is the side test, now angle-robust).
- The largest-component fallback (θ→90° exactly) is gated at 1e-6; production
  hardening with sorted-global-vertex-ID orientation (option c) remains optional.

---

## Frontera Dc2 re-run (job 7746601, commit 800b328) — side fix NECESSARY but NOT sufficient; dip-channel blow-up remains

CORRECTION (user, 2026-05-22): the blow-up is **NOT fixed**.  It would still
occur at **t≈1.0 s**; the run aborts EARLY at **t≈0.455 s** only because the
hard R-101 DOFData-consistency guard trips first.  The side-determination fix
(removing the 71.6° side-flip) is necessary and real — the run now advances
cleanly to t≈0.455 s with `V_max=0` and R-101 passing at machine precision
(1.28e-15 → 1.75e-14) — but the underlying **dip-channel (V1/tau1) cross-rank
inconsistency + structural instability is still the driver of the t≈1.0 s
blow-up**.  The guard converts a t≈1.0 s numerical blow-up into a clean t≈0.455 s
abort on the first detectable inconsistency.  **Resolution was never the issue**;
the side-flip was one cause, the dip channel is the other.  Decision (user): the
dip channel MUST be made cross-rank-consistent AND numerically stable — the
strike-slip-only `SEAS_FORCE_V1_ZERO` lock is NOT acceptable (SAFS has genuine
dip slip).

It then trips a DIFFERENT, pre-existing failure the R-101 guard catches:

```
R-101 ... FAILED.  Field 'V1' at centroid (607518, 3706359, -4543)
differs by 1.6368 across the two ranks (rel_diff=1.0).  [t≈0.455 s, nucleation patch]
```

### Diagnosis: the DIP channel (V1/tau1) is exercised for the FIRST time
- `SolveLSW_TPV205` (`tpv205_friction.hpp:146-165`): `V1 = V_abs·tau1_total/|τ|`,
  `tau1_corr = tau1_trial − η_s·V1`.  So **V1=0 exactly ⟺ tau1_total=0 exactly**.
  One rank has dip total-traction ≈0 (V1=0), the other nonzero (V1=1.637).
- The R-101 verifier reports the MAX-rel-diff field (`wave_operator.inl:5893`),
  so `tau1_corr` ALSO diverged (rel≈0.29, ~7.6 MPa = η_s·1.637) but V1's rel=1.0
  is larger, hence V1 is named.  Both are dip-channel quantities → the divergence
  is in `tau1_total = tau1_0 + tau1_nuc + tau1_trial` (dip), not strike.
- **Why never seen before:** ALL benchmarks (TPV102/104/205, BP5) are pure
  strike-slip → V1≡0 by construction.  `fault_face_flux.cpp:211-218` documents a
  KNOWN structural instability in this very channel ("closed-loop V1 amplification
  …1-ULP perturbation grows ~1.07/step…correct only for strike-slip") whose only
  mitigation is the strike-slip-only env-var lock `SEAS_FORCE_V1_ZERO`.  SAFS is
  the first config with GENUINE dip slip, so it cannot use the lock and is the
  first to actually run the dip channel.

### Leading hypotheses (need a local dip-slip reproducer to distinguish)
- **H-A (most likely): dip-prestress `tau1_0` projected on the per-rank STORED
  `tangent1`** (which is sign-flipped across a shared face: `t1_A = −t1_B`,
  v61 design) **instead of the canonical `can_t1`** → `tau1_0` inconsistent
  across ranks → `tau1_total` diverges.  Clean, fixable (project prestress onto
  the canonical frame).
- **H-B: the documented structural V1 instability** amplifies an unavoidable
  1-ULP cross-rank seed (different FP op-order on the two ranks) by ~1.07/step.
  Formulation-level; harder.
- **H-C: `can_t1` not actually bit-identical across ranks.**  Less likely (the
  tilted test exercised the canonical frame — but only with STRIKE nucleation, so
  the dip channel was never validated cross-rank).

### Proposed next step
Extend `test_rupture_tilted_fault_serial_vs_parallel` (the tool that nailed the
side bug) to drive DIP slip — nucleate in `tau1_nuc` and/or add a dip-shear
prestress — and assert cross-rank `V1` consistency at np=2, printing per-rank
`tau1_0 / can_t1 / tau1_total`.  That distinguishes H-A/H-B/H-C locally before
any change to the friction code ("Files Requiring Extreme Care").
