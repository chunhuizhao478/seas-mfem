# TPV31 benchmark deviation — findings & improvement report

Date: 2026-06-02
Author: Claude (code-debug)
Runs analysed (comparison PNGs, MFEM-vs-SCEC overlay):
  * `tpv31/plot_out_p1_aderO2_normal_7761406`  (p1 · ADER-O2 · matrix Riemann)
  * `tpv31/plot_out_p2_aderO3_normal_7761407`  (p2 · ADER-O3 · matrix Riemann)
Reference: `tpv31/benchmark_data/scec_{eqdyna,seisol}/` (30 on-fault stations each)
Spec: `tpv31/benchmark_document/TPV31_32_Description_v03.pdf`

> **No code was changed.** This is a diagnosis + recommendations report
> (per CLAUDE.md: report findings and ask before applying fixes).

> **UPDATE 2026-06-02 (TPV102 control + concrete root cause).** The same σ_n
> collapse + spurious dip-slip appears in the **TPV102 spatial run**
> (`tpv102/plot_out_p1_aderO2_pu_normal_7761452`, job 7761452) — which is
> **homogeneous + scalar Godunov + pure upwind + rate-state**, i.e. it shares
> *none* of TPV31's material/flux/law features. The only common factor is the
> **`seas_spatial_dyn_driver` at 500 ranks**. Tracing it pinned the concrete
> cause: **the spatial driver does NOT apply the fault-locality partition that
> the native `tpv102/104/205_driver.cpp` all use** (`fault_locality_partition.hpp`,
> `--partition-fault-locality`). Without it ParMETIS cuts the y=0 fault plane
> across ranks, and the +y/−y per-sub-step ghost exchange seeds the
> `[[v_n]]`/σ_n leak. See the new §5.1.

---

## 1. Bottom line

**The on-fault normal stress `σ_n` (station column 8 = `DOFData::sigma_n_corr`)
drifts strongly tensile over the run — from its correct depth-proportional
initial value to a near-uniform ~22 MPa — when it must remain *exactly
constant*. This is a numerical instability in the elastodynamic/wave-operator
layer, NOT a problem-setup or friction-law error.** It is the **same
"speckle/sliver blow-up" the SAFS spatial-driver investigation already diagnosed**
(`safs/.../spatial_dynamic_rupture_sliver_blowup_2026-05-25.md`): a **seed** —
a spurious, secularly-growing normal-velocity jump `[[v_n]]` from the per-sub-step
ADER predictor / shared-face ghost exchange (R-1303/R-1601), seeded at mesh
slivers — amplified by the **LSW `max(σ_n,0)` free-slide** into a fault-opening
runaway. (My first pass blamed a rank-2 stress *rotation*; the May investigation
**refuted** that — `project_safs_vn_leak_not_frame`.) Every other deviation
(shear undershoot, over-slip, spurious dip slip, secondary pulses, and the
~11 s truncation) is downstream of this single cause. TPV31 is hit hard because
it is **LSW** (the `max(σ_n,0)` amplifier; rate-state TPV104 only leaked ~0.5 MPa),
runs at **500 ranks** (many shared/seam fault faces), past the **reflection
window** with **no PML**, and its config enables **none** of the existing
mitigations (`sigma_n_strength_floor`, `SEAS_NOOPENING`).

The problem **setup is spec-exact and reference-confirmed** (see §4).

---

## 2. Observed deviations (symptom)

Comparing MFEM (solid) against EQdyna (dashed) and SeisSol (dotted):

| # | Quantity | SCEC reference | MFEM (p1 & p2, identical) | Verdict |
|---|----------|----------------|---------------------------|---------|
| D1 | **σ_n at 7.5 km** | **constant −60.6305 MPa for all 15 s** (literally identical at t=0.005 and t=15.0 in both EQdyna & SeisSol) | starts 60.6 MPa, **collapses to ~22 MPa** beginning ~t≈5–6 s | **non-physical — bug** |
| D2 | σ_n vs depth | constant at each depth (24.5 / 60.6 / 81.1 MPa at 0 / 7.5 / 10 km) | starts correct, all drift toward ~22 MPa | bug (same cause) |
| D3 | strike shear τ at 7.5 km | settles ~25.7 MPa | undershoots to ~22 MPa | consequence of D1 |
| D4 | down-dip slip | ~5×10⁻⁶ m (zero) | grows to mm-scale, p1≠p2 | consequence of D1 |
| D5 | strike slip (final) at 7.5 km | 2.17 m | ~2.3–2.5 m (over-slip) | consequence of D1 |
| D6 | V_strike secondary pulses | single clean pulse | extra bumps at t≈6–10 s | consequence of D1 |
| D7 | run duration | full 15 s | **traces stop ~11 s** | consequence of D1 (see §3) |

The σ_n collapse is **roughly simultaneous across all stations** (onset ~5–6 s,
worsening to ~11 s) — it is a global secular drift, not a propagating physical
phase. Both p1/ADER-O2 and p2/ADER-O3 show it (so it is not a pure
discretization-order effect); the p1↔p2 divergence on the dip channel (D4)
indicates the drift is seeded/grown numerically.

---

## 3. Why σ_n *must* be constant, and the runaway it triggers

TPV31's fault is the planar `y=0` plane. The 1-D velocity/stress structure
varies **only with depth** → the material is **identical on both sides of the
fault** at every depth. For a pure strike-slip rupture on such a fault there is
no mechanism to change the fault-normal stress; `σ_n(t) ≡ σ_n0(depth)`. The
SCEC reference traces confirm this to the digit (column 8 is bit-constant).

The fault-flux code itself encodes this: `FaultFaceFlux::EvaluateADER_LSW`
asserts equal impedance across the fault
(`fault_face_flux.cpp:778`, `homog_ok(Zp_plus,Zp_minus)`), and for TPV31 it
holds — so the **fault Riemann solve is the *same code path validated by
TPV205***. The "bimaterial" character of TPV31 lives entirely in the **volume**
(depth-layered medium), handled by `BimaterialWaveOperator`.

**Runaway feedback (why TPV31 explodes while TPV104 only leaked ~0.5 MPa):**
the LSW strength is `τ_str = μ_eff·max(σ_n,0) + C₀`
(`tpv205_friction.hpp:174–177`, default `sigma_n_floor = 0`). A spurious
tensile drift in σ_n directly lowers `τ_str`, which lets the fault slip more,
which radiates more, which feeds the leak — a positive feedback. As
σ_n → 0 the strength → C₀ (= 0 below 2.4 km), the fault free-slides, slip rate
spikes, and the run blows up / stalls. **This explains D7** (the ~11 s
truncation; whether it is a NaN-abort, CFL blow-up, or wall/step limit should
be confirmed from the Frontera `*.out/*.err` for jobs 7761406/7761407 — these
logs are **not** synced into the local plot dirs). TPV104 uses *regularized
rate-and-state*, which has no hard `max(σ_n,0)` threshold, so its identical
leak stayed bounded at ~0.5 MPa.

---

## 4. Problem-setup audit vs spec — PASS

Every physics parameter in `tpv31/configs/tpv31.toml` (and `_p2`) was checked
against the spec PDF and cross-checked against the reference traces:

| Item | Spec | Config | Status |
|------|------|--------|--------|
| Geometry | 30 km × 15 km vertical fault, reaches surface, hypocenter (0,7.5 km) | y=0 plane x∈[-15,15] km z∈[-15,0] km, hypo (0,0,-7500) | ✓ |
| Stress tensor | σ₁₁=σ₃₃=−60·(μ/μ₀), σ₁₃=+30·(μ/μ₀), rest 0, μ₀=32.03812032 GPa | rotated → σ_xx=σ_yy=60, σ_xy=30 (compression+) | ✓ (rotation verified component-by-component) |
| **σ_n0 / τ0 at 7.5 km** | 60.63 / 30.32 MPa | matches; **EQdyna/SeisSol report σ_n=−60.6305 MPa → exact match** | ✓ |
| Material 1-D profile | discontinuities at 2400/5000/10000 m, graded 2400→5000 | 5-layer `depth_profile_1d`, R-007 tie-break | ✓ |
| Nucleation | +4.95·(μ/μ₀) MPa, r≤1400, cosine taper to r=2000, pure right-lateral | `instantaneous_overstress_circular` radius 1400, taper 600, Δτ 4.95e6, μ-scaled | ✓ |
| Friction | μ_s=0.58, μ_d=0.45, d_0=0.18 m | matches | ✓ |
| Cohesion | 0.000425 MPa/m·max(0,2400−depth) (1.02 MPa surface→0 at 2.4 km) | depth-taper rule, grad 425 Pa/m, ref_depth 2400 | ✓ |
| Resolution | 50 m | 50 m on fault; cohesive zone Λ₀≈650 m ≈ 13 elems @ 50 m → adequate | ✓ |
| Final time | 15 s | tfinal=15 s, dt_max=0.01 s | ✓ |

**Setup is not the cause of the deviations.** Two *minor* setup weaknesses
worth noting (they do not explain a 38 MPa monotonic drift but are secondary
contributors and should be closed once D1 is fixed):

* **Absorbing boundary is basic upwind-with-zero-background**
  (`godunov_flux.cpp:455`, `AbsorbingTotal` → `Interior` with `Q_bg=0`), not a
  Clayton–Engquist characteristic absorber or PML. It absorbs normal incidence
  well but **leaks for oblique P/S**. PML exists but is **off** (`use_pml=false`).
  The box (50 km) clears the worst-case reflection-free distance (48.75 km) by
  only **2.6 %**; the nearest leaky reflections (x-walls / bottom, 35 km from the
  fault) return at **t≈10.8 s** — close to where the trace dies, so reflected
  energy is plausibly a *late* aggravator, but the σ_n onset at ~5–6 s is too
  early to be explained by it.
* **ADER hot-path optimizations are always on** (`--deriv-cache
  --shared-ck-recursion`). Their equivalence is verified "at round-off" by a
  unit test that (per the job header) is the homogeneous-material equivalence
  test — **not** validated on the heterogeneous matrix path. A control run
  *without* these flags is a cheap way to rule them out.

---

## 5. Root-cause localization (where the σ_n leak lives)

Ruled **out** (with evidence):
* Friction / fault Riemann solve — shared with validated TPV205; `σ_n_corr =
  σ_n0 + σ_n_nuc + σ_n_trial` with `σ_n_trial` from the standard Godunov
  formula (`fault_face_flux.cpp:62`). ✗ not the cause.
* Per-DOF impedance init — set once, both sides equal (`spatial_setup.hpp:85`),
  constant in time. ✗.
* Stress / nucleation / material setup — spec-exact, reference-confirmed (§4). ✗.
* Interior bimaterial volume flux — uses correct per-element material on each
  side, no `(1,1,1)` placeholder leak (`bimaterial_wave_operator.inl:365–421`,
  `godunov_flux_bimaterial.cpp:71–296`). ✗.
* "Slip-pinned-at-fault-border" not implemented — *not a bug*: the fault is a
  **finite surface patch** (`Physical Surface(101)` only on x∈[-15,15],
  z∈[-15,0]); slip tapers to zero at the border geometrically because the
  surrounding faces are welded interior faces, not frictional. The comment at
  `tpv31_stations.hpp:64` refers to this. ✗.

Ruled **in** — the root cause is a **two-layer seed + amplifier**, already
diagnosed (and partially fixed) by the SAFS spatial-driver investigation
(same driver as TPV31). The authoritative analysis is
`safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_sliver_blowup_2026-05-25.md`
and `…/document/DEBUG_speckle_normal_velocity_jump_2026-05-23.md`
(root cause confirmed 2026-05-24, with a passing local np=2 test).

> **Correction to my first pass:** the rank-2 stress-**rotation/frame** theory
> (carried over from the *April* TPV104 doc) was **explicitly refuted** by the
> later, more thorough May investigation —
> `project_safs_vn_leak_**not_frame**`: the frame canonicalization
> (`sign_flipped`, `wave_operator.inl:2365`) is clean; `[FRAME]` orthonormality
> is ~16 ULP by construction. The leak is **not** the rotation.

**Layer 1 — SEED: a spurious normal-velocity jump `[[v_n]]` from the per-sub-step
ADER predictor / ghost exchange on SHARED (MPI partition-seam) fault faces**
(R-1303/R-1601). The friction iterator reads the element-local Cauchy–Kovalevskaya
**point predictor** `Q̃(τ_o)` at each sub-step; on a shared face the neighbour
side comes from a per-sub-step **ghost exchange** (`wave_operator.inl:2030–2174`,
the `qpd.sign_flipped` gather at :2362). The predictor's normal-velocity jump
**diverges from the time-integrated macro solve by ~10⁸** (`DEBUG_speckle…`:
predictor `dv_n` ≈ 3.8e6 m/s vs macro trial jump 0.015 m/s at the seed), and the
divergence is **secular** across macro steps (gain ~1.3/step) — a positive-
eigenvalue spatial mode, **not** a CFL/time-stepping problem (`sliver_blowup`
§6b: dt is already worst-sliver-limited by the insphere-CFL). At the seed
`[[v_n]] ≈ 5.2 m/s ≈ 41 % of the tangential slip rate`, resolution-independent.
**Mesh slivers set *where*** (worst-conditioned faces seed first; for this
geometry, the shallow z=0 fault-surface trace and the bimaterial depth
interfaces are the suspects), **accumulated reflection energy sets *when*** (no
PML; box reflections return ~10.8 s, surface breakout ~4–5 s).

**Layer 2 — AMPLIFIER: the LSW `τ_str = μ·max(σ_n,0) + C₀` free-slide**
(`tpv205_friction.hpp:174`), field-confirmed in `sliver_blowup` §2c. The growing
`[[v_n]]` drives σ_n tensile; the moment `σ_n_total ≤ 0` the strength collapses
to `C₀` (= 0 below 2.4 km), the Godunov flux opens the fault, shear traction →0,
and slip is unbounded (`V = |τ|/η_s`). More slip → more `[[v_n]]` → more tensile
σ_n → **positive-feedback runaway** (σ_n "chatter": tensile-open ↔
compressive-at-μ_d each macro-step; the 0.1 s station/ParaView cadence
**aliases** the per-step flips into the apparent smooth "collapse to ~22 MPa").
This is why a bounded ~0.5 MPa leak (TPV104, rate-state, *no* `max(σ_n,0)`
threshold → no runaway) becomes a total collapse + ~11 s blow-up in TPV31 (LSW).

**Status of fixes / why TPV31 is still hit:**
* `550c640` landed an **R-004 fix** (correct byNODES ghost-exchange neighbour
  read on shared fault faces) — addresses one seed, but the SAFS `sliver_blowup`
  (05-25, *after* it) still blew up, so the instability is not fully closed.
* `sigma_n_strength_floor` **is implemented** (`fault_face_flux.hpp:513`,
  applied at `tpv205_substep_iterator.cpp:126`, `friction_substep_iterator.cpp:184`,
  `fault_face_flux.cpp:233`) and the `SEAS_NOOPENING` cap exists
  (`tpv205_friction.hpp:153`) — **both default-disabled, and TPV31's config sets
  NEITHER**, so the Layer-2 amplifier runs unchecked.
* TPV31 **already uses `mixed_flux = "none"`** (pure upwind, dissipative), so the
  SAFS "central-flux removes dissipation" destabilizer (`sliver_blowup` §6b) does
  **not** apply here — TPV31's residual seed must be the shared-face predictor/ghost
  `[[v_n]]` path and/or sliver/bimaterial-interface modes under upwind.

> Confirming *which* residual seed dominates for TPV31 needs the existing
> `SEAS_DIAG_SLIP` / `[FRAME]` / `[MACRO]` / `SEAS_DIAG_BLOWUP` traces on a short
> Frontera slice (local full-mesh runs are disallowed). The mechanism class is
> confirmed; the per-config dominant seed is not yet pinned for TPV31.

### 5.1 TPV102 control experiment + the concrete root cause (2026-06-02)

The TPV102 spatial run (`tpv102/plot_out_p1_aderO2_pu_normal_7761452`, job 7761452)
is the decisive control. It exhibits the **identical signature** as TPV31:
σ_n collapses from 120 MPa to **~20 MPa** at the hypocenter (DRDG3D reference is
flat at −120.0000→−120.005 MPa over 15 s), spurious dip slip grows to **0.32 m**
(reference ~8×10⁻⁷ m), and shear undershoots. Yet TPV102 is:

| feature | TPV31 | **TPV102 (control)** | ⇒ rules out |
|---|---|---|---|
| material | bi-material (depth 1-D) | **homogeneous** | the matrix/bimaterial path |
| interior flux | matrix Riemann | **scalar Godunov** | the bimaterial flux |
| flux type | pure upwind | **pure upwind** | central flux (already off in both) |
| friction law | **LSW** | **rate-state (aging)** | the LSW law / rotation theory |
| driver / ranks | spatial / 500 | **spatial / 500** | — (the ONLY shared factor) |

So the seed is **independent of material, interior-flux kind, and friction law** —
it is purely the **`seas_spatial_dyn_driver` shared-fault path at np>1**.

**Onset tracks the local rupture front** (not a global clock or boundary
reflection): at the TPV102 hypocenter σ_n departs at ~1.5 s; at the strike-9 km
station it stays flat at 120 MPa until ~5 s (rupture arrival) then collapses to
~40 MPa by 12 s. The collapse depth scales with **time-since-local-rupture**
(hypocenter: ~10 s → 20 MPa; 9 km: ~7 s → 40 MPa) — the signature of a **secular
leak generated by the actively-slipping fault** on shared faces.

**Why the difference in outcome (amplifier is law-dependent):**
- **RS (TPV102/104):** strength uses `|σ_n|` (`fault_face_flux.cpp:233`), so a
  reduced σ_n only weakens the fault — no free-slide. The σ_n leak collapses σ_n
  by ~100 MPa (83 %) and drives mm-to-dm spurious dip slip, but the run stays
  **bounded** and completes 12 s (V_strike ~4 m/s). **Wrong, not divergent.**
- **LSW (TPV31):** `max(σ_n,0)` (`tpv205_friction.hpp:174`) → once the leak drives
  σ_n tensile, strength → 0, free-slide, **runaway → blow-up ~11 s.**

This is the definitive answer to the SAFS open question ("would rate-state blow
up at the same slivers? — the σ_n/`[[v_n]]` leak is law-agnostic"): **the leak is
law-agnostic and full-magnitude in RS; only the *runaway* is LSW-specific.**

**THE CONCRETE ROOT CAUSE — missing fault-locality partition in the spatial driver.**
`dynamic/fault_locality_partition.hpp` documents this exact failure (written for a
TPV104 np≥4 dip-pollution event, debug doc 2026-04-25_pm §13):

> "ParMETIS subdivides along Y, **cutting the y=0 fault plane**; ghost-cell
> exchanges between +y and −y ranks introduce a topology-induced FP
> non-associativity that the rupture amplifies … to mm-scale slip_dip."

The fix (`BuildFaultLocalityPartitioning` / `VerifyFaultLocality`) forces both
elements of every fault face onto the same rank, so the fault plane is never cut
→ no +y/−y seam → no seed. **All three native drivers apply it**
(`tpv104_driver.cpp:39,485,1184–1201`; `tpv102_driver.cpp`; `tpv205_driver.cpp`
all `#include "../dynamic/fault_locality_partition.hpp"` and offer
`--partition-fault-locality`). **`spatial_dyn_driver.cpp` does NOT** — it has no
include, no flag, and no call (grep: the only "partition" hit, line 1178, is an
unrelated comment). So every spatial-driver dynamic run (TPV31, TPV102/104/205-
spatial, SAFS) at np≥4 cuts the fault plane and carries this seed; the native
drivers do not. **This is the bug, and it is config-/driver-level, not physics.**

---

## 6. Recommended next steps (in priority order)

**A. THE FIX — port fault-locality partitioning into the spatial driver (removes the seed).**
 1. Wire `dynamic/fault_locality_partition.hpp` into `drivers/spatial_dyn_driver.cpp`
    exactly as the native drivers do (`tpv104_driver.cpp:39,485,1184–1201`):
    build the serial mesh, collect fault faces (attr 101), call
    `BuildFaultLocalityPartitioning(...)`, hand the partition vector to the
    `ParMesh` ctor, and `VerifyFaultLocality(...)` (assert 0 violations). Add a
    `--partition-fault-locality` flag (default **on** for the dynamic spatial
    driver — there is no reason to cut a fault plane). This keeps both elements
    of every fault face on one rank → no +y/−y seam → **no `[[v_n]]` seed**. It
    fixes **accuracy** (σ_n stays constant; dip slip → 0) for **all** spatial-
    driver dynamic runs (TPV31, TPV102/104/205-spatial, SAFS) — not just the
    blow-up. This is the real fix; everything below is confirmation or a band-aid.

**B. Confirm decisively (cheap):**
 2. Re-run TPV102-spatial (or a short TPV31 slice) at 500 ranks **with**
    fault-locality partitioning. Prediction: σ_n stays flat at 120 / 60.6 MPa,
    dip slip → ~0, TPV31 reaches 15 s. (If a quick check is wanted before the
    port: run at **np=1** — no shared faces — and confirm both are clean.)
 3. The existing `[FRAME]`/`[MACRO]`/`SEAS_DIAG_SLIP`/`SEAS_DIAG_BLOWUP` traces
    (commits `3926d7b`, `97dd189`) on a ~1 s slice quantify the residual seed if
    any remains after the partition fix.
 4. **Pull the Frontera `*.out/*.err`** for 7761406/7761407 (TPV31) to record the
    ~11 s termination cause, and confirm the binary contained `550c640`.

**C. Band-aids (do NOT fix accuracy — use only to get a completed TPV31 run):**
 * `sigma_n_strength_floor_pa = 10e6` under `[friction]` (`fault_face_flux.hpp:513`)
   caps `τ_str` and **stops the LSW blow-up**, but σ_n still leaks — it merely
   converts TPV31 into TPV102's "bounded-but-wrong" state. Useless for RS
   (TPV102/104 use `|σ_n|`, already bounded). **Not an accuracy fix.**
 * Freezing σ_n to its per-DOF background (generalize `SEAS_TPV104_FREEZE_SIGMA_N`,
   `fault_face_flux.cpp:126`) masks the symptom but leaves the dip-slip / `[[v_n]]`
   pollution. Prefer (A).

**D. Secondary (independent):** PML / larger box (`use_pml=false`, 2.6 % margin)
   to remove reflection energy that *aggravates* the seed once it exists.

**C. Secondary (independent of the leak fix):**
 * Turn on PML (`use_pml=true`) or shrink the box + add a PML sponge, OR switch
   the absorbing flux to a characteristic absorber, to remove the marginal
   boundary-reflection contribution and the 2.6 % box margin.

**D. Add a regression guard (general improvement):**
 * For planar-fault, symmetric-material configs (TPV31/205/102/104), assert/warn
   when `max_t |σ_n_corr − σ_n0|` on the fault exceeds a small tolerance
   (e.g. a few % of σ_n0). σ_n is analytically constant for this whole problem
   class, so this is a cheap, high-signal tripwire that would have caught this
   immediately. Wire it into the station writer or a per-macro-step reduction.

---

## 7. Evidence index

* σ_n constant in reference: `benchmark_data/scec_{eqdyna,seisol}/tpv31_*_x2_0_x3_7.5.txt`
  col 8 = −60.6305 MPa at both t=0.005 s and t=15.0 s.
* σ_n collapse in MFEM: `plot_out_p2_aderO3_normal_7761407/tpv31_faultst000dp075.png`
  (and dp000, dp100, st120dp075) — bottom-left panel.
* Trial-traction formula: `dynamic/fault_face_flux.cpp:49–85`.
* LSW strength `μ·max(σ_n,0)+C₀`: `dynamic/tpv205_friction.hpp:166–182`.
* Equal-impedance fault assertion: `dynamic/fault_face_flux.cpp:771–784`.
* Per-QP rotation: `dynamic/wave_operator.inl:2072–2098`;
  `dynamic/godunov_flux.cpp:304` (`BuildFrame`).
* Prior precedent + recommended fix: `debug_document/tpv104_debug_document/tpv104_sigma_n_leak_root_cause_2026-04-25.md`.
* Existing mitigations/diagnostics: commits `156ea34`, `97dd189`, `3926d7b`.
