# TPV31 — residual on-fault DIP-SLIP excess (after σ_n is controlled)

Date: 2026-06-06
Author: Claude (code-debug)
Runs analysed (MFEM-vs-SeisSol overlay PNGs):
  * `tpv31/plot_out_p2_aderO3_pu_50m_hybrid_7771023`  (p2 · ADER-O3 · pure upwind · matrix Riemann · **hybrid y-symmetric mesh**)
  * `tpv31/plot_out_p1_rk45_mf_normal_7774534`        (p1 · RK45 · **mixed/central flux** · matrix Riemann · unstructured mesh)
Reference: `tpv31/benchmark_data/scec_seisol/` (SeisSol, 30 on-fault stations)
Git: `b09dbfc` (branch `system/spatial_dyn_driver`)

> **No code was changed.** Diagnosis + recommendation only (per CLAUDE.md:
> report findings and ask before applying fixes). Verification requires a
> Frontera A/B (no local full-mesh runs).

> **Relationship to the prior doc** (`tpv31_sigma_n_drift_findings_2026-06-02.md`):
> that doc diagnosed the **catastrophic σ_n collapse → LSW blow-up** (jobs
> 776xxxx, σ_n 60→22 MPa, run dies ~11 s). That failure mode is **GONE** in
> these 777xxxx runs: σ_n is flat (only ±0.005–0.01 MPa chatter) and both runs
> reach 15 s. The prior doc listed dip slip (its "D4") as *downstream of the
> σ_n collapse*. **This doc updates that:** with σ_n now controlled, a small
> residual dip-slip excess **survives** — so it is NOT merely downstream of σ_n.
> It is the **dip-tangential component of the same antisymmetric-mode leak**,
> with its own seed that the σ_n fixes do not close.

---

## 1. Bottom line

The on-fault **dip slip** (`slip1`, down-dip component) is biased away from
the SeisSol reference. It is **tiny in magnitude** (~1×10⁻⁴ m at the surface
vs ~2 m strike slip — i.e. ~5×10⁻⁵ of the strike slip), so this is an
**accuracy-floor / polish** issue, not a correctness-threatening bug. Strike
slip, strike stress, normal stress, and the dip slip *at depth* all match well.

**Root cause:** a spurious **dip-tangential velocity jump `[[v_dip]] = [[v_t1]]`**
across the fault, driven by **3-D mesh asymmetry** and **amplified at the free
surface**, which the dissipative (irreversible) friction **rectifies** into a
secular dip-slip offset. It is the exact tangential analog of the on-fault
σ_n drift — same trial-traction equation, different component (see §4).

> **§9 (2026-06-06) added a LOCAL UNIT-TEST PROOF of this mechanism** (extended
> `test_fault_planar_serial`) AND a partial cure test. Headline: the y-mirror
> mesh annihilates the σ_n/`[[v_n]]` leak (~10¹²×) but leaves the `[[v_dip]]`
> leak intact (ratio ~0.9) — the TPV31 paradox reproduced at unit scale. The
> cure question is **not** settled: **over-integration alone does NOT reduce the
> dip leak (it slightly increases it)**, and the resample half of the dealiasing
> is untestable in this harness. See §9 — and treat §7's over-integration
> recommendation as **superseded by §9**.

It does **not** reduce with polynomial order because it is an
**aliasing/asymmetry** artifact, not a resolution error (§5). Its sign and
magnitude are **method-dependent** (upwind-p2-hybrid: +1.2×10⁻⁴ m at surface;
mixed-flux-p1: −3.5×10⁻⁵ m) because the two runs close *different* sub-channels
of the same leak (§6).

---

## 2. Symptom (from the comparison PNGs)

3rd-row-left panel = "Slip Dip (m)"; black = MFEM, light red = SeisSol.

| station | depth | MFEM dip slip (upwind p2 hybrid) | MFEM dip slip (mixed-flux p1) | SeisSol ref |
|---|---|---|---|---|
| dp000 | 0 km (free surface) | **+1.2×10⁻⁴ m** (step at rupture front, holds) | **−3.5×10⁻⁵ m** | ~0 |
| dp010 | 1 km | **−1.0×10⁻⁴ m** (sign flips vs surface) | — | ~0 |
| dp075 | 7.5 km | +2×10⁻⁵ m | (small) | ~0.5×10⁻⁵ m |
| dp120 | 12 km | +3×10⁻⁵ m | — | ~2×10⁻⁵ m (real) |

Key signatures (all load-bearing):
1. **Surface-concentrated.** Largest at z=0, decays into the bulk. At depth
   (≥7.5 km) MFEM ≈ reference (the small *physical* dip slip from the depth-
   varying medium).
2. **Depth-oscillating near the surface.** +1.2×10⁻⁴ at z=0, −1.0×10⁻⁴ at
   z=−1 km, back to +2×10⁻⁵ at depth. A free-surface standing-/Rayleigh-wave
   fingerprint, **not** a constant offset.
3. **Ratchets.** Steps up at local rupture-front arrival and holds flat — the
   hallmark of a dissipative friction integrating a small dip traction into a
   permanent offset (no reversal).
4. **Method-dependent sign/magnitude** (upwind-hybrid + vs mixed-flux −).
5. **p-independent** (p1 vs p2 on the SAME hybrid mesh,
   `tpv31_p{1,2}_pureupwind_50m_hybrid`, both show it).
6. μ_eff weakens cleanly 0.58→0.45 → the friction law itself is healthy.

---

## 3. Why this is NOT a frame/sign coding bug (ruled out)

A constant coordinate-frame or sign error in the dip channel would be
**uniform in depth and identical across methods**. Instead the artifact is
surface-concentrated, depth-oscillating, and flips sign between upwind and
mixed-flux runs, while MFEM and the reference **agree** on the small physical
dip slip at depth (both ≈ +2×10⁻⁵ m). The dip-channel code is the standard
exact-Riemann formula and is correct (§4). Ruled out.

---

## 4. Mechanism — the dip leak is the σ_n leak in another component

`FaultFaceFlux::ComputeTrialTraction` (`dynamic/fault_face_flux.cpp:49–85`)
computes all three fault channels identically from the **velocity jump**
(fault-local frame: t1 = dip, t2 = strike; Eq. 7a–c):

```
sigma_n_trial = eta_p * ( [[v_n]]    + sigma_n^+/Zp^+ + sigma_n^-/Zp^- )   (7a)  -> normal
tau1_trial    = eta_s * ( [[v_dip]]  + tau1^+/Zs^+    + tau1^-/Zs^-   )    (7b)  -> DIP
tau2_trial    = eta_s * ( [[v_str]]  + tau2^+/Zs^+    + tau2^-/Zs^-   )    (7c)  -> strike
```
(`[[a]] = a^- - a^+`.) The LSW solve (`dynamic/tpv205_friction.hpp:188–205`)
then sets the slip rate **exactly parallel to the total shear traction**:

```
V_abs = max(tau_abs - tau_strength, 0)/eta_s
V1(dip)    = V_abs * tau1_total / tau_abs        slip_dip   += V1 * dt   (irreversible)
V2(strike) = V_abs * tau2_total / tau_abs
```

So a spurious dip-tangential velocity jump `[[v_dip]]` → nonzero `tau1_trial`
→ nonzero `V1` → **accumulating dip slip**. This is mathematically the **same
equation** as the σ_n channel (7a vs 7b), so whatever asymmetry leaks into
`[[v_n]]` (the documented σ_n drift) also leaks into `[[v_dip]]`. They are one
leak in two components.

**Why it peaks at the free surface:** the dip-tangent t1 ≈ the vertical (z)
direction on this vertical fault. The free surface permits large vertical
particle motion (Rayleigh/SV reflection), so the z-velocity field — and hence
any +y/−y sampling asymmetry in `[[v_dip]]` — is **largest near z=0**. Below,
the medium constrains v_z and the leak decays to the physical value.

This matches the Phase-0 harness result
(`tests/unit/test_fault_planar_serial.cpp`): the σ_n drift is a *rectified
accumulation of a sustained velocity-jump leak* that **vanishes to round-off on
a mirror-symmetric mesh and is ~270 kPa on an asymmetric mesh** (ratio ~2.6×10¹²)
— mesh asymmetry is NECESSARY. The dip slip is that same leak in the
dip-tangential component.

---

## 5. Why raising polynomial order does NOT help

The leak has two asymmetry-gated sub-channels, **both aliasing-class**:
- **Friction aliasing.** The nonlinear LSW strength `μ(δ)·σ_n + C₀` and the
  slip-parallel decomposition are evaluated **pointwise at the fault QPs**
  (degree `2*order`, `--fault-overint 0`). Under-integrating a nonlinear flux
  aliases high modes into the resolved traction — including a dip-tangential
  component on an asymmetric mesh. Raising `p` adds **more** high modes to
  alias; it does not dealias. (This is exactly the SeisSol "resample"/over-
  integration story — the cure is over-integration, not higher order.)
- **Near-fault upwind dissipation** (Zhang 2023): the upwind Riemann flux is
  dissipative; on an asymmetric mesh it dissipates the +y/−y halves unequally,
  seeding an antisymmetric jump. Also geometry-driven, p-independent.

Both are properties of the **mesh geometry + quadrature**, not the polynomial
resolution. Hence p1≈p2 on the same mesh, as observed.

---

## 6. Why the two runs differ in sign/magnitude (and why each σ_n fix leaves dip slip)

Both runs tame σ_n (the `[[v_n]]` / **y-antisymmetric** channel) by different
means, but neither closes the dip (`[[v_dip]]` / **z-tangential**) channel:

| run | σ_n fix | what it closes | why dip slip survives |
|---|---|---|---|
| upwind p2 **hybrid mesh** (7771023) | y-mirror prism strip near fault | the **y-antisymmetric** `[[v_n]]` mode → σ_n flat | the strip is **y-mirror only, not z-symmetric**; the z-tangential `[[v_dip]]` is driven by the unstructured/3-D + free-surface asymmetry the strip does not remove |
| **mixed/central flux** p1 (7774534) | central flux removes upwind dissipation | the **upwind-dissipation** sub-channel | central flux does **not** dealias the **friction-aliasing** sub-channel → residual `[[v_dip]]` |

The opposite signs are then expected: different residual seeds (y-mesh-residual
vs friction-aliasing) on different meshes.

This also yields a **decisive A/B** (§7): the mixed-flux run has *already*
removed upwind dissipation, so if **over-integration** then removes its
residual dip slip, friction-aliasing is confirmed as the dip seed; if dip slip
persists, the residual is the **volume/mesh-asymmetry wave-field** coupling
(point to near-surface mesh quality / a z-structured near-fault strip).

---

## 7. Recommended verification (Frontera A/B — ask before launching)

All knobs already exist; nothing here is a code change.

1. **Over-integration on the mixed-flux run (the cleanest test).** Re-run
   `tpv31_p1_rk45_mixedflux_50m_normal` with `--fault-overint 2` (raises the
   fault-face quadrature to degree `2*(order+2)`; dealiases the friction·flux
   solve). Central flux already removed upwind dissipation, so this isolates
   the friction-aliasing seed. **Prediction:** surface dip slip shrinks
   markedly if friction aliasing is the residual seed.
   - Note: for **LSW**, `--fault-resample` resamples slip **magnitude only**
     (direction preserved — `friction_substep_iterator.cpp:284–314`), so it
     will *not* by itself remove a directional dip artifact. Use
     `--fault-overint`, with `--fault-resample` optional/secondary.
2. **Over-integration on the upwind hybrid run.** Add `--fault-overint 2` to
   `tpv31_p2_pureupwind_50m_hybrid`. Tests whether dealiasing also reduces the
   z-residual on the y-symmetric mesh.
3. **Diagnostic:** if available, dump `max_F |[[v_dip]]|` and the dip-channel
   `tau1_trial` near z=0 on a short slice — confirm `[[v_dip]]` tracks the
   dip-slip ratchet (the analog of the Phase-0 `[[v_n]]` probe). Consider
   extending `test_fault_planar_serial` to also report `max_F|[[v_t1]]|` so the
   dip leak has a serial, symmetric-vs-asymmetric unit-level testbed (it should
   vanish on the mirror mesh, like `[[v_n]]`).

Secondary (independent, from the prior doc, §C/§D): `use_pml=false` and the
2.6 % absorbing-box margin add oblique-reflection energy that can *aggravate*
the dip seed once it exists; the dip onset tracks the rupture front (not the
~10.8 s reflection return), so reflections are at most a late aggravator, not
the seed.

---

## 8. Evidence index

* Trial traction (3 channels from velocity jump): `dynamic/fault_face_flux.cpp:49–85`.
* LSW strength + slip-parallel decomposition: `dynamic/tpv205_friction.hpp:166–205`.
* Over-integration wiring (`--fault-overint`, off by default):
  `drivers/spatial_dyn_driver.cpp:563–567, 1317–1346`.
* Resample (LSW = magnitude-only): `drivers/spatial_dyn_driver.cpp:1448–1515`;
  `dynamic/friction_substep_iterator.cpp:284–314`.
* Phase-0 mesh-asymmetry result: `tests/unit/test_fault_planar_serial.cpp`
  (symmetric mesh → `[[v_n]]`/σ_n ≈ round-off; asymmetric → ~270 kPa).
* Job scripts: `jobs/tpv31_spatial/upwind_ader/tpv31_p2_pureupwind_50m_hybrid_normal.sbatch`
  (FAULT_OVERINT=0, hybrid mesh), `jobs/tpv31_spatial/mixedflux_rk45/tpv31_p1_rk45_mixedflux_50m_normal.sbatch`
  (`--mixed-flux adjacent`, unstructured mesh, no over-integration).
* Prior σ_n doc: `tpv31_sigma_n_drift_findings_2026-06-02.md`.

---

## 9. LOCAL UNIT-TEST PROOF (2026-06-06) — mechanism confirmed, over-integration-alone refuted

`tests/unit/test_fault_planar_serial.cpp` was extended to measure, alongside the
existing `max_F|[[v_n]]|` (normal) probe, the **dip-tangential jump
`max_F|[[v_dip]]|`** (= `[[v_z]]` on the y-normal slab) and the strike jump
`[[v_strike]]`. The slab is a homogeneous strike-slip nucleation pulse, so the
ideal dip slip is **zero** → any `[[v_dip]]` is pure leak. The "symmetric"
variant is y-mirror-symmetric near the fault (it does NOT remove the z-direction
Kuhn-dicing asymmetry — exactly like TPV31's y-mirror hybrid mesh). Serial p1,
3 over-integration levels. **Runs locally (~3 min each); no Frontera, no full
mesh.**

### Proof matrix (peak over the run)

| `--fault-overint` | dip `[[v_dip]]` asym | dip sym | **dip ratio asym/sym** | normal `dSigN` asym | normal sym | **normal ratio** |
|---|---|---|---|---|---|---|
| 0 (baseline) | 7.39e-3 m/s | 8.80e-3 | **0.84** | 2.70e5 Pa | 1.19e-7 | **2.27e12** |
| 1 | 8.96e-3 | 9.54e-3 | 0.94 | 3.18e5 | 1.34e-7 | 2.37e12 |
| 2 | 9.59e-3 | 10.25e-3 | 0.94 | 3.01e5 | 1.04e-7 | 2.89e12 |

(For scale, the real strike rupture `[[v_strike]]` ≈ 7.0–7.5e-2 m/s, so the dip
leak is ~10–13 % of the genuine slip rate — small but not round-off.)

### What it PROVES (decisive)

1. **The dip leak is a distinct, y-symmetry-IMMUNE channel.** The y-mirror mesh
   drives the normal/σ_n leak to round-off (ratio ~2–3×10¹²) but barely touches
   the dip leak (ratio ~0.84–0.94) — at **every** over-integration level. This
   is the unit-scale reproduction of the TPV31 paradox: the y-symmetric hybrid
   mesh flattened σ_n yet left the dip slip, because dip slip lives in the
   z-tangential channel the y-mirror does not constrain. **Confirms §4/§6.**
2. **Mesh asymmetry is the seed** (the dip leak exists with NO free surface in
   this harness — the free surface only AMPLIFIES it in TPV31).

### What it REFUTES / leaves open (be honest)

3. **Over-integration ALONE is NOT the dip-slip cure.** Raising the fault-face
   quadrature (k=0→2) monotonically **increases** the dip leak (+30 %) and does
   not reduce the normal leak either. So the earlier §7 "turn on `--fault-overint`"
   recommendation is **not supported** by the controllable local test.
4. **BUT the full SeisSol dealiasing (over-integration + RESAMPLE) is NOT
   testable in this harness:** `--fault-resample` is inert here, and the harness
   calls `wave.AdvanceADER` directly (it does **not** route through
   `RateStateSubStepIterator`/`LinearSlipWeakeningIterator` where the production
   resample lives). Over-integration is only the *enabling* half (per
   `[[seisol-resample-is-degree-N-not-Nminus1]]`: "resample is identity without
   over-integration; over-int first, resample second"). So the full cure remains
   **unproven**, not disproven.
5. The harness regime (weak decaying rate-state pulse, peakV ~0.12 m/s)
   **under-represents friction-aliasing-of-a-vigorous-front**, which TPV31's LSW
   rupture has. So the harness leans toward the **near-fault upwind-dissipation**
   channel (Zhang 2023, geometry-driven) as the dominant dip seed *here* — a
   channel over-integration cannot fix but **central/mixed flux can** (and the
   mixed-flux TPV31 run did show a smaller surface dip slip, 3.5e-5 < 1.2e-4).

### Corrected recommendation (supersedes §7)

- **Do NOT bank on `--fault-overint` for dip slip.** Local proof shows it alone
  worsens the leak. Also note `SetFaultOverint(k>0)` is **guarded incompatible
  with the mixed-flux path** (`wave_operator.inl:634–639`, Phase-1 scope), so it
  cannot even be applied to the mixed-flux production run.
- **Two real candidate cures, ranked by current evidence:**
  1. **Central/mixed flux near the fault** removes the upwind-dissipation seed
     the harness implicates; the mixed-flux TPV31 run already has the smaller
     surface dip slip. Lowest-risk lever available today.
  2. **Full dealiasing (over-int + resample)** on the **pure-upwind** path —
     still a candidate for the friction-aliasing component of a vigorous front,
     but it must be tested with BOTH knobs, either by (a) a Frontera A/B on
     `tpv31_p2_pureupwind_50m_hybrid` with `--fault-overint 2 --fault-resample`,
     or (b) extending the harness to route through the substep iterator so the
     resample is live locally (the honest local test of the full cure).
- **Mesh quality** in the dip-amplification (near-surface) zone is a third lever
  but a z-symmetric mesh is impractical with a free surface at z=0.

### Harness changes (committed to the test, not the solver)

`test_fault_planar_serial.cpp`: generalized `MaxVnJumpInterior` →
`MaxVCompJumpInterior(comp)`; added per-step `[[v_dip]]` peak tracking +
`[[v_strike]]` reporting; added two regression guards codifying the finding
(dip leak present on asym; y-mirror does NOT suppress it, unlike σ_n). All
prior acceptance checks still pass (now 8 tests). No solver code changed.
</content>
</invoke>
