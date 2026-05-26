# Spatial dynamic rupture — late sliver-seeded V_max blow-up (job 7748818, 2026-05-25)

Third documented SAFS dynamic-rupture blow-up. **Slip-weakening (LSW)**, not
rate-and-state. Distinct in *timing and location* from the two prior events but
shares the **mesh-sliver** root cause with the 2026-05-20 north-tip event.

| | this event (2026-05-25) | 2026-05-22 speckle | 2026-05-20 north-tip |
|---|---|---|---|
| job | 7748818 | 7744292 | (postevent) |
| config | `..._slip_weakening_safs_projected_stress_resolution_Dc2.toml` | `..._slip_weakening_safs_projected_stress.toml` | — |
| onset | **late, t≈28 s** (runaway t≈42 s) | early, t≈0.84 s (nucleation) | post-event, t≈3 s |
| where | **mid-fault x≈365 km + shallow trace x≈535 km** | on/next to nucleation patch | far north tip |
| recovers? | **NO** → 5e64, never recovers | yes → back to ~3 m/s | no → NaN |
| root cause | **mesh slivers** (this doc) | self-limiting on-patch speckle | mesh slivers (shallow trace) |

> **Correction to the 2026-05-22 doc.** That doc states the `500m z0embed`
> mesh is "sliver-free." **It is not.** This investigation found the
> single worst triangle in the entire fault mesh (min-angle **9.1°**, tri_q
> **0.245**) at x≈365 km, plus a sliver band at the shallow surface trace.
> The 05-22 *speckle* was indeed not sliver-driven (it was on the
> well-meshed patch and recovered); but the mesh as a whole still carries
> sliver patches, and **this** blow-up sits squarely on them.

Run header (from `/Users/chunhuizhao/Documents/spatial_dyn_resDc2_7748818.log`):

```
config:  safs/.../spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml
mesh:    safs/.../meshing/results/msh/safs_fault_box_nwcut_500m_lcfar3000_z0embed.msh
law:     slip_weakening        ader order: 2     cfl: 0.5      use pml: no
tfinal:  100 s                 ranks: 500        num_fault_global = 128151
[time]   dt = 0.000349884 s    nsteps = 285809
[derived] mesh h_min = 197.029 m (scalar cp)
[derived] L_nuc/h_min min/max = 6.09 / 19.45   (NB: should be >= 10)   <-- under-resolved warned
[derived] PASS: initial conditions are well-posed.
```

---

## 1. Symptom — late, super-exponential, non-recovering

The rupture nucleates and propagates **correctly** for ~25 s, then a localized
instability grows and runs away:

| t (s) | step | V_max (m/s) | regime |
|---|---|---|---|
| 0–1.75 | 0–5k | 0 → 4 | nucleation ramp (gradual_overstress) |
| ~3 | 8.4k | ~10 | rupture breaks out, propagates |
| 3 → 28 | 8k–80k | **stable 10–12** | healthy front; `n_rupturing` 3.6k→8.8k; **survives the 6.94 s reflection window** |
| 28 | 80k | 15 | first departure from the plateau |
| 33 / 42 | 94k/120k | 30 / 100 | slow climb |
| 44 / 46 | 125k/132k | 1e3 / 1e6 | **super-exponential cascade** |
| 64.6 | 184.6k | **5e64** | fully blown up; `n_rupturing`=113k/128k (fault-wide) |

The LSW physics, fault flux, and nucleation are therefore **correct** — a clean
~10 m/s rupture for ~25 s. The blow-up is a **late, growing numerical
instability**, not a setup error.

---

## 2. Trace — normal-stress "chattering" (fault open/close), not physical slip

Tracing the dominant runaway DOF `(365100, 3.83892e6, z=−2488 m)` through the
cascade (`[DIAG-ONSET]` records, `SEAS_DIAG_BLOWUP=1`):

```
V=749     sigma_n_corr=-7.15e8   tau2_corr=-2.3e-7    <- tensile sigma_n, traction ~ 0 (fault OPEN)
V=1027    sigma_n_corr=+1.39e8   tau2_corr=-6.6e7     <- compressive again, traction back
V=19896   sigma_n_corr=-2.6e10   tau2_corr=-7.9e-6
V=53044   sigma_n_corr=+9.2e10   ...
```

`sigma_n_corr` **flips sign every step (tensile↔compressive) with exponentially
growing amplitude**; whenever it goes tensile, the Godunov fault flux opens the
fault and the shear traction collapses to ~0 → no restoring friction → V runs
away. This is a grid-scale **normal-stress oscillation**, the same `[[v_n]]`
class noted in `project_safs_vn_leak_not_frame` (per-sub-step predictor/ghost,
R-1303/R-1601) — here amplified to destruction.

---

## 2c. Field-data confirmation (fault.vtkhdf, 2026-05-26) — tensile σ_n + free slip CONFIRMED

The full fault output (`/Users/chunhuizhao/Documents/fault.vtkhdf`, 10.4 GB, 1301 snapshots
t∈[0, 64.59 s], 128 151 points) lets us verify the σ_n-chatter / fault-opening mechanism from
the *saved fields* (`normal_stress`, `traction_{strike,dip}`, `slip_rate_{strike,dip}`),
independent of the log. Figure: `blowup_clusterA_sigma_n_tau_V_2026-05-25.png`.

Cluster A point (idx 114460, x365029 y3838938 **z=−2560 m**) trajectory:

| t [s] | σ_n [Pa] | \|τ\| [Pa] | \|V\| [m/s] | τ/\|σ_n\| | state |
|---|---|---|---|---|---|
| 24–43 | +8.42e7 | 2.56e7 | **0** | 0.304 | locked at the μ_d threshold |
| 44.5 | +5.69e7 | 3.29e7 | 1.88 | 0.58 | instability arrives; τ overshoots |
| 45.0 | **+1.02e9** | 6.3e7 | 275 | 0.06 | σ_n spikes 20× (oscillation explodes) |
| 45.4 | **−2.91e9** | 1.77e9 | 1744 | 0.61 | **first TENSILE — fault opens** |
| 45.6 | −3.23e10 | 1.35e9 | 1.8e4 | **0.04** | tensile → τ collapsing (free slip) |
| 46.0 | −4.70e11 | 1.96e11 | 1.25e5 | 0.42 | super-exponential cascade |
| 64.0 | **−1.27e69** | 1.82e66 | 5.76e62 | **0.001** | tensile, τ→0 (free slip) |

**Both hypotheses confirmed:**
- **Tensile normal stress — YES.** σ_n at the blow-up points goes strongly **negative** (compression
  is +; first tensile at t≈45.3 s; 199/1301 snapshots tensile), reaching **−1.3e69 Pa**. The
  healthy reference point (nucleation patch, idx 62662) stays **compressive the whole run**
  (σ_n ∈ [49.5, 59.5] MPa, 0/1301 tensile, |V|max 3.6 m/s) — the contrast is decisive.
- **Free slip (τ→0) — YES, on the tensile half-cycles.** When σ_n is tensile the fault opens and
  τ/|σ_n| collapses to **0.001–0.04** (vs μ_d=0.30) — it sheds shear restoring traction. It is a
  **normal-stress *chatter*** (σ_n alternates tensile-open ↔ compressive-at-μ_d each macro-step;
  the 0.05 s snapshots alias the per-step flips the log showed), so the median τ/|σ_n| over the
  whole runaway (0.33–0.39) is the mixture of free-slip (~0) and at-limit (~0.30) half-cycles.

**Sequence (mechanism).** A locked patch (V=0, τ/σ_n=0.30) → the slowly-growing instability
arrives (~t=44.5 s) → the σ_n oscillation amplitude explodes over ~5 macro-steps (84 MPa→1 GPa)
→ σ_n crosses **tensile** (~t=45.4 s) → fault **opens**, τ→0 (free slip) → no restoring force →
super-exponential |V| runaway. The Godunov fault flux *correctly* opens the fault under tension
(faults can't sustain tension) — but here the tension is **numerical** (the spurious oscillation),
so the opening converts a normal-stress wiggle into a real loss of shear restoring → positive
feedback. **This is the runaway stage; the σ_n oscillation that grows it (t=28–45 s) is the
spatial-operator instability of §6b (central-flux suspect).**

**Spatial spread** (snapshot counts, of 128 151 points): tensile σ_n / |V|>50 m/s grow from
**801 / 0** at t=40 s → 1179 / 39 at t=46 s → 12 780 / 14 682 at t=52 s → **48 343 / 111 574**
at t=64.6 s. So it seeds at ~800 sliver/near-surface points and the opening+radiation
destabilizes neighbours until it is fault-wide — matching the log's `n_rupturing → 113 k`.

## 3. Where — two clusters, both on low-quality (sliver) fault elements

`[DIAG-ONSET]` logs the **global argmax DOF per step** (printed only when
V_max>10). Aggregating peak |V| and argmax-hit-count per location gives 9
runaway hotspots (peak |V| > 50 m/s; physical peak was ~10–15), in **two tight
clusters, both far from the nucleation patch** (x≈607 km, z=−5.1 km):

| cluster | x | y | z | peak \|V\| (m/s) | argmax hits |
|---|---|---|---|---|---|
| **A** (mid-fault) | 365 040–365 100 | 3.8389e6 | **−2.5 to −2.7 km** | up to **8.9e64** | 31339 + 17258 + … |
| **B** (shallow trace) | 535 056–535 364 | 3.7553e6 | **−25 to −122 m** | up to **7980** | 8799 + 2527 + … |

(Full list: `…/meshing/results/vtu/safs_fault_blowup_points.csv`.)

---

## 4. ROOT CAUSE — the blow-up elements ARE mesh slivers

`meshing/code/locate_fault_slivers.py` + a fault-wide quality sweep on
`safs_fault_box_nwcut_500m_lcfar3000_z0embed.msh` (fault physical tag 101,
42 358 triangles):

| metric | fault-wide median | fault-wide worst | the 20 marked blow-up cells |
|---|---|---|---|
| min angle | 58.9° | **9.10°** | median 32.8°, min **14.0°** |
| tri_q (1 = equilateral) | 0.9995 | 0.245 | median 0.795, min **0.312** |
| min edge | 499 m | 112 m | down to **112 m** |
| aspect | 1.02 | 4.43 | up to 3.86 |

Correlation:
- The marked cells' **median** min-angle (32.8°) is at the **0.80th percentile**
  fault-wide — worse than 99.2 % of all fault triangles.
- The **single worst triangle in the whole mesh** (9.10°, tri_q 0.245,
  centroid `(365012, 3838896, −4277)`) is at **cluster A**; the 2nd and 3rd
  worst (10.6°, 11.8°) are also at cluster A.
- Fault-wide only 274/42358 (0.65 %) triangles have min-angle < 30°, yet
  **8 of the 20** marked cells do — a ~40× enrichment.
- Cluster A: min-angle 9–14°, min_edge **112 m**. Cluster B (shallow trace):
  min-angle 20–25°, min_edge **181 m**. The well-meshed fault is ~equilateral
  with 500 m edges.

### Mechanism (where × when × how)
The fault is excellent almost everywhere, but sliver patches have an effective
element size (min_edge 112–181 m, sliver **inradius** far smaller) well below
the `h_min = 197 m` used to set the global explicit `dt`. Those slivers are the
worst-conditioned points on the fault. Sequence:
1. **Where** = slivers. Cluster A is the global-worst triangle; cluster B is the
   shallow-trace sliver band. Both are exactly where the runaway localizes.
2. **When** = energy accumulation. The run is stable through the first reflection
   passes (t 7–28 s); with **no PML** (only 1st-order absorbing BCs on attrs
   103/104) and `tfinal = 100 s ≫ 6.94 s` reflection time, reflected-wave energy
   accumulates and at t≈28 s the marginally-unstable sliver modes begin to grow.
3. **How** = σ_n chatter. The growing mode drives `sigma_n` tensile → fault opens
   → shear traction → 0 → unbounded slip-rate growth (§2).

Slivers set *where*; accumulated reflection energy + the σ_n/`[[v_n]]` leak set
*when/how*. Cluster A being the literal global-worst sliver, and the dominant
runaway, is not coincidence.

This matches `locate_fault_slivers.py`'s own provenance note (written for the
2026-05-20 event): "spurious north-tip velocity coincided with sliver triangles
at the shallow fault-surface trace."

---

## 5. Artifacts produced (2026-05-25)

In `…/meshing/results/vtu/`:
- `safs_fault_blowup_marked.vtu` — fault surface (42 358 tris) with cell arrays
  `is_blowup` (0/1), `blowup_peakV`, `blowup_loghits`, **`min_angle_deg`,
  `tri_q`, `aspect`, `min_edge_m`**. Color by `min_angle_deg` / `tri_q`: the
  marked blow-up cells sit in the dark sliver patches.
- `safs_fault_blowup_points.vtu` / `.csv` — the 982 argmax hotspots; `is_blowup`
  flags the 9 with peak |V| > 50 m/s.

Generators (re-runnable, `conda activate pythonenv`):
- `…/meshing/code/mark_blowup_on_fault.py` — log → marked VTU/points (`VBLOW=50`).
- `…/meshing/code/locate_fault_slivers.py` — spatial quality probe
  (`--target X Y Z --radius R --fault-tag 101`).

Reproduce the cluster-A quality:
```
python3 code/locate_fault_slivers.py \
  results/msh/safs_fault_box_nwcut_500m_lcfar3000_z0embed.msh \
  --target 365100 3838920 -2487.95 --radius 1500
```

---

## 6. Recommended actions (highest leverage first)

1. **Repair the sliver patches** at x≈365 km (z −0.9 to −4.3 km) and the shallow
   surface trace (x≈535 km, z≈−0.1 km). This is the *where*; remove it and the
   instability has no seed. Plans already exist:
   `meshing/docs/PLAN_mesh_quality.md`, `PLAN_retriangulate.md`,
   `EXPLORE_tmop_mesh_optimization.md`. The "z0embed" surface-embedding step is
   the likely sliver source (fault meeting the z=0 free surface at a grazing
   angle → thin triangles).
2. **Confirm the local-CFL story on the 3D tets** (TODO). The surface triangle is
   a proxy; the explicit `dt` is set by the *tet* inradius/cp. Check the tet
   quality at cluster A — if a sliver tet's stable dt < 3.5e-4 s, that is a true
   local-CFL violation. (Open: why stable for 28 s then runaway? — accumulated
   energy, not an immediate hard CFL trip.)
3. **Run inside the valid window** regardless: `tfinal ≤ ~7 s` or wire real PML
   (the driver's `use_pml` flag is currently a no-op — see REVIEW.md R-026). The
   blow-up is entirely in the already-reflection-contaminated regime (t>6.94 s).
4. **Survey other-code sliver handling** (next task): how SeisSol treats
   degenerate/sliver elements (local time stepping, element dropping, quality
   gates, CFL per-element). Inform whether MFEM-SEAS needs per-element dt or a
   mesh-quality hard gate.

## 6b. UPDATE (2026-05-25, later) — it is NOT a time-stepping/CFL problem

Verified `ComputeMaxDt` (wave_operator.inl:5634) and the tet geometry:

- MFEM-SEAS sets a **single global dt = min over tets of `cfl_eff · h_e / cp_e`**,
  where `h_e` = the tet **inscribed diameter** (`6V/A_total`, wave_operator.inl:97)
  — the *same* geometric measure SeisSol uses (`GlobalTimestep.cpp:computeCellTimestep`,
  insphere). So MFEM-SEAS is **not** naive; the dt is already sliver-aware.
- Tet inscribed-diameter over the 3.69 M tets: **min 42.0 m**, p1 127 m, median 322 m.
  The reported `h_min = 197 m` is a *different* metric (`pmesh.GetElementSize(e,0)`,
  used only for the L_nuc/h display) — the **dt actually uses 42 m**:
  `dt = (cfl/9)·0.9·42/cp ≈ 3.5e-4 s` matches the log exactly (cfl 0.5, DG `3(2p+1)`
  factor, mixed-flux 0.9).
- The cluster-A blow-up **tets** have min inscribed-diameter **46 m** = 1.10× the global
  limiter (42 m). Cluster B = 95 m (2.27×). **So dt is sized essentially for the very
  elements that blow up** — they are CFL-stable by construction.

**Conclusion: dt is already set by the worst sliver, and they still blow up ⇒ this is a
SPATIAL-operator instability, not a temporal/CFL one.** A CFL-stable dt only guarantees
stability if the semi-discrete operator is dissipative/neutral (all eigenvalues ≤ 0). If the
DG fault/interior operator has a weakly **positive** eigenvalue on the distorted sliver
elements, the mode grows as exp(λt) for **any** dt — no dt reduction and no LTS fixes it.
The "stable 28 s then exp runaway at a CFL-safe dt" signature is exactly a positive-eigenvalue
spatial mode fed by accumulating energy, not a CFL trip.

**Leading suspect — non-dissipative central flux.** This run used `mixed_flux = adjacent`:
the log shows **200 308 interior faces on central (non-dissipative) flux** (Zhang 2023). The
`ComputeMaxDt` comment itself flags "central flux is non-dissipative." Central flux + a 9°
sliver's poor conditioning removes the dissipation that would otherwise damp the σ_n/`[[v_n]]`
leak (`project_safs_vn_leak_not_frame`, R-1303/R-1601 predictor/ghost). SeisSol uses **upwind
(Godunov)** interior flux — dissipative — so its insphere-CFL dt is sufficient; MFEM-SEAS's
central-flux choice removes that safety net on distorted elements.

**Reflections set the *when*, not the *root*.** No PML + tfinal 100 s ≫ 6.94 s lets reflected
energy accumulate; the unstable mode needs amplitude to become visible (onset t≈28 s). PML /
shorter tfinal would *delay* it, not cure a positive-eigenvalue mode.

### Answers to the two design questions
**Q1 — sliver standard.** The project gate already exists (`check_mesh_quality.py`: Q1 min tet
edge ≥ 100 m, Q2 Joe-Liu eta > 0.1). The 500 m mesh **passes** it (emin 112 m, eta_min 0.184)
yet blows up — so the gate is **too lenient for this scheme**. Define the standard on the
quantity that drives each failure mode: (a) **fault triangle** min-angle ≥ ~25–30° / tri_q ≥
~0.5–0.6 (the σ_n leak lives on the fault face; the blow-up tris were 9–25°); (b) **tet**
Joe-Liu eta floor raised to ~0.2–0.3 (0.1838 is marginal). But note: tightening the mesh gate
*reduces the amplifier*; it does not, by itself, make a non-dissipative operator stable.
**Q2 — time stepping.** Not the lever. dt is already worst-sliver-limited; reducing it or
adding SeisSol-style LTS (which is a *performance* optimization, same insphere-CFL per cluster)
will **not** stop a positive-eigenvalue spatial mode. The fixes are spatial: remove the slivers
(mesh repair), and/or restore dissipation (`mixed_flux = none` / pure upwind), and/or penalize
the fault-normal `[[v_n]]` weld.

### Discriminating experiment (recommended, from a checkpoint near t≈27 s, step ~78000)
1. Rerun with **`mixed_flux = none`** (pure upwind). If the blow-up vanishes or moves vastly
   later → central flux is confirmed as the instability source. (Cheapest, highest-information.)
2. Rerun with **PML / tfinal = 6 s**. If onset just delays/weakens → reflections are an
   amplifier, not the root.
3. Rerun on a **sliver-repaired** mesh (TMOP / retriangulate the x≈365k + shallow-trace patches)
   at the same flux. If stable → the sliver conditioning was the amplifier.

## 7. Open questions
- Is `h_min = 197 m` the bulk-tet min, or does it miss the 112 m sliver edges?
  (The fault min_edge 112 < 197, so `h_min` does not see the smallest fault
  feature — worth auditing how `dt_cfl` is computed.)
- Would the rate-and-state path blow up at the same slivers? (RS regularizes
  friction but the σ_n/`[[v_n]]` leak is law-agnostic.)
- Does the `--print-derived` gate need a **mesh-quality** check (min-angle /
  per-element CFL) in addition to `L_nuc/h_min`?

---

# Implementation Plan: optional compressive normal-stress strength floor (σ_n cap on shear strength)

> **STATUS: DRAFT — for /code-review → /code-fix → /code-implement.** Author 2026-05-26.
> Targets the confirmed free-slip half of the blow-up (§2c): under tension the LSW strength
> collapses to 0 (`SolveLSW_TPV205` clamps σ_n to 0), so the fault free-slides. This adds an
> **opt-in** floor so the strength stops being proportional to σ_n below a threshold.

## Overview
Expose a **compressive normal-stress floor** `σ_n_floor` (TOML, default = disabled). In every
shear-strength computation, replace the normal stress used for strength with
`σ_n_strength = max(σ_n_total, σ_n_floor)`. Above the floor the strength is the usual
`μ·σ_n` (proportional); **below the floor it is the constant `μ·σ_n_floor`** (no longer
proportional to σ_n, and never the spurious tensile `0` (LSW) or `|σ_n|` (RS)). Set
`σ_n_floor = 10 MPa` in the SAFS configs. This breaks the σ_n→strength→radiation feedback that
drives the runaway (§2c, §6b) while staying **byte-exact** for the TPV/BP5 regressions when
disabled.

**Governing change (math).** For both laws the friction strength is `τ_strength = μ_·σ_n^{str}`
where the friction coefficient `μ_` is unchanged (LSW: `μ_eff(δ)`; RS: `f_V = a·asinh(V·C)`), and
```
σ_n^{str} = max(σ_n_total, σ_n_floor)        (floor enabled, σ_n_floor ≥ 0)
σ_n^{str} = max(σ_n_total, 0)   [LSW] / |σ_n_total| [RS]   (floor disabled — current behavior)
```
Only the normal stress entering the **strength** is floored; the written-back `σ_n_corr`
(the `normal_stress` output channel) is untouched. `σ_n_floor = 10 MPa` ⇒ below 10 MPa
compression the strength saturates at `μ·10 MPa` (a cohesion-like constant), so a spurious
tensile excursion can no longer drive free slip.

## Constraints
- **Byte-exact regression is non-negotiable.** TPV205/102/104 and BP5 must be bit-identical when
  the floor is disabled. ⇒ the option is **opt-in via a disabled sentinel** (`σ_n_floor < 0`),
  and when disabled each path keeps its *exact* current expression (LSW `max(σ_n,0)`, RS
  `|σ_n|`). The benchmarks set no floor ⇒ untouched.
- **`SolveLSW_TPV205` is the byte-exact oracle kernel** (tpv205_friction.hpp). The new parameter
  MUST default to the current behavior (`σ_n_floor = 0.0` ⇒ `max(σ_n, 0.0)`), so the standalone
  TPV205 driver + the oracle iterator that don't pass it are unchanged.
- **CLAUDE.md "don't change sign conventions / don't revert a fix."** This is *additive*: the
  existing `max(σ_n,0)` (LSW), `abs(σ_n)` (RS), and the `SEAS_NOOPENING` env guard all stay; the
  floor is a new gated branch. Cite §2c/§6b as justification in the commit.
- **Single carrier.** All four strength sites are reached through one `FaultFaceFlux` (the
  iterators bind `FaultFaceFlux&`), so the floor lives as one `FaultFaceFlux` member, set once
  from config; no per-DOF plumbing in v1 (flag per-DOF as a follow-up).
- **Floor is a compressive (positive) stress**; `μ_` is the *active* coefficient (continuous at
  the threshold). See "Decisions for /code-review" for the `μ_s`-literal alternative.

## Strength sites to change (all four; mapped 2026-05-26)
| # | file:line | path | current σ_n in strength |
|---|---|---|---|
| 1 | `dynamic/tpv205_substep_iterator.cpp:110` | **SAFS LSW** (LswFrictionIterator→Tpv205) — the blow-up path | via `SolveLSW_TPV205` → `max(σ_n,0)` |
| 2 | `dynamic/fault_face_flux.cpp:804` | LSW `EvaluateADER_LSW` (shared/inline) | via `SolveLSW_TPV205` → `max(σ_n,0)` |
| 3 | `dynamic/fault_face_flux.cpp:921` | LSW forced-rupture variant | via `SolveLSW_TPV205` → `max(σ_n,0)` |
| 4 | `dynamic/fault_face_flux.cpp:226` and `:562` | RS substep / shared (aging) | `abs(σ_n)` |

## Phase 1 — config plumbing (no physics change)

### Goal
A `σ_n_floor` value parses from TOML into the config and reaches `FaultFaceFlux`; default
disabled (sentinel `< 0`); all existing tests byte-exact.

### Files to Modify
- `spatial/code/spatial_friction.hpp` — add to `struct SpatialFrictionConfig`:
  `real_t sigma_n_strength_floor_pa = -1.0;  ///< <0 = disabled; >=0 = floor [Pa] on σ_n in shear strength`.
- `spatial/code/spatial_friction.cpp` — in the `[friction]`-level parse (the function that fills
  `SpatialFrictionConfig`, alongside `law`), read optional key
  `sigma_n_strength_floor_pa` (default -1.0). Validate: if present, `>= 0` and finite, else MFEM_ABORT.
- `dynamic/fault_face_flux.hpp` — add private member `real_t sigma_n_strength_floor_ = -1.0;`
  and public setter `void SetSigmaNStrengthFloor(real_t v) { sigma_n_strength_floor_ = v; }`.
- `drivers/spatial_dyn_driver.cpp` — after `FaultFaceFlux fault_flux(...)` (`:1258`):
  `fault_flux.SetSigmaNStrengthFloor(cfg.sigma_n_strength_floor_pa);`.

### Interfaces
- `SpatialFrictionConfig::sigma_n_strength_floor_pa` (real_t, default -1.0).
- `FaultFaceFlux::SetSigmaNStrengthFloor(real_t)` / member `sigma_n_strength_floor_`.

### Edge Cases
- key absent → -1.0 (disabled). key present and `< 0` → MFEM_ABORT (a negative floor is the
  sentinel, not a value). key `= 0` → allowed (≡ LSW `max(σ_n,0)`; for RS switches `abs`→`max(·,0)`).

### Acceptance Criteria
- [ ] New unit `test_spatial_friction_config`: a TOML with `sigma_n_strength_floor_pa = 10.0e6`
      parses to `10e6`; absent → `-1.0`; negative value aborts.
- [ ] `make test` byte-exact (no behavior consumes the field yet).
- [ ] `seas_spatial_dyn_driver` builds; `--dry-run` prints the floor (add to the banner).

### Dependencies
Depends on: nothing. Required by: Phases 2, 3.

## Phase 2 — apply the floor in the LSW strength (the SAFS blow-up path)

### Goal
LSW shear strength uses `max(σ_n, σ_n_floor)`; default-disabled is bit-identical to today;
enabled (10 MPa) caps the strength so a tensile σ_n no longer free-slides.

### Files to Modify
- `dynamic/tpv205_friction.hpp` — `SolveLSW_TPV205(...)`: add a **trailing defaulted** parameter
  `real_t sigma_n_floor = 0.0` (after `tau2_corr`). Replace line 157
  `const real_t sigma_n_pos = std::max<real_t>(sigma_n_total, 0.0);`
  with `const real_t sigma_n_pos = std::max<real_t>(sigma_n_total, sigma_n_floor);`.
  Default `0.0` ⇒ `max(σ_n,0)` ⇒ byte-exact for every caller that omits it.
- `dynamic/fault_face_flux.cpp:804` and `:921` — pass `sigma_n_floor_arg` (see below) as the new
  trailing arg to `SolveLSW_TPV205`.
- `dynamic/tpv205_substep_iterator.cpp:110` — pass `flux_.SigmaNStrengthFloorForLSW()` as the new
  trailing arg (the iterator holds `FaultFaceFlux& flux_`).
- `dynamic/fault_face_flux.hpp` — add accessor
  `real_t SigmaNStrengthFloorForLSW() const { return sigma_n_strength_floor_ >= 0.0 ? sigma_n_strength_floor_ : 0.0; }`
  (maps the disabled sentinel `-1` → `0.0`, i.e. current LSW `max(σ_n,0)`).

### Interfaces
- `SolveLSW_TPV205(..., real_t &tau1_corr, real_t &tau2_corr, real_t sigma_n_floor = 0.0)`.
- `FaultFaceFlux::SigmaNStrengthFloorForLSW() const`.

### Edge Cases
- floor disabled (`-1`) → accessor returns `0.0` → `max(σ_n,0)` → byte-exact.
- floor `= 10 MPa`, σ_n tensile (−2.9 GPa) → `sigma_n_pos = 10 MPa`, `τ_strength = μ_eff·10 MPa`
  (finite) → `V_abs = (|τ| − μ_eff·10 MPa)/η_s` is bounded (vs free slip `|τ|/η_s`).
- the `mu_eff >= 0.5·mu_s_barrier` lock (line 117) and `SEAS_NOOPENING` (line 144) branches are
  **before** line 157 and are unchanged.

### Acceptance Criteria
- [ ] New unit `test_lsw_strength_floor` on `SolveLSW_TPV205` directly: with
      `σ_n=−2.9e9, μ_eff=0.3, η_s=4.6e6, |τ|=1.7e9`: floor `0.0` → `τ_strength=0`,
      `V_abs=|τ|/η_s` (free slip); floor `10e6` → `τ_strength=0.3·10e6=3e6`,
      `V_abs=(|τ|−3e6)/η_s` (bounded, ≈ free-slip minus a small finite cap) — assert the floored
      `V_abs < free-slip V_abs` and `τ_strength == μ_eff·σ_n_floor`.
- [ ] `seas_test_*` for TPV205 / the substep-iterator parity **byte-exact** with floor disabled
      (the default-`0.0` arg path).
- [ ] `make test` green.

### Dependencies
Depends on: Phase 1. Required by: Phase 4.

## Phase 3 — apply the floor in the RS strength (aging path), gated

### Goal
RS shear strength uses `max(σ_n, σ_n_floor)` when enabled; disabled keeps `abs(σ_n)` (byte-exact
TPV102/104).

### Files to Modify
- `dynamic/fault_face_flux.cpp:226` and `:562` — replace `std::abs(s.sigma_n_total)` /
  `std::abs(sigma_n_fric)` with a helper:
  ```cpp
  const real_t sn_str = (sigma_n_strength_floor_ >= 0.0)
                        ? std::max(s.sigma_n_total, sigma_n_strength_floor_)   // (sigma_n_fric at :562)
                        : std::abs(s.sigma_n_total);
  real_t strength = sn_str * f_V;
  ```
  Note the friction *solver* call (`solver_.Solve(..., std::abs(s.sigma_n_total), ...)` at `:201`/`:548`)
  also takes σ_n — **decide in /code-review** whether to floor that argument too (it sets the Brent
  bracket / `Theta`); v1 floors only the explicit `strength` to minimize the byte-exact surface.

### Edge Cases
- disabled (`-1`) → `abs(σ_n)` (current). enabled → `max(σ_n, floor)`.
- TPV102/104 σ_n is compressive ~120 MPa ⇒ with no floor set, `abs ≡` current ⇒ byte-exact.

### Acceptance Criteria
- [ ] `seas_test_tpv102_nuc_callback_parity`, `test_resolve_rate_state_guards`, the RS np=2 test:
      **byte-exact** with floor disabled.
- [ ] Unit: RS strength at `σ_n=−1e9` with floor `10e6` → `strength = f_V·10e6` (bounded),
      not `f_V·1e9`.

### Dependencies
Depends on: Phase 1. Required by: Phase 4.

## Phase 4 — enable on the SAFS configs + blow-up acceptance

### Files to Modify
- `safs/.../config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml`
  (and `_Dc8/_Dc10`, and the RS configs) — add under `[friction]`:
  `sigma_n_strength_floor_pa = 10.0e6   # cohesion-like floor: below 10 MPa compression, strength = mu*10MPa (constant)`.

### Acceptance Criteria
- [ ] `--dry-run` banner reports `sigma_n strength floor: 10 MPa`.
- [ ] **Blow-up arrest (cluster, the real test):** rerun the LSW Dc2 case (same as job 7748818,
      `mixed_flux=adjacent`) with the floor at 10 MPa to ~t=35 s; assert V_max stays bounded
      (e.g. < 1e3 m/s) past the original t=28 s onset — i.e. the free-slip runaway is broken.
      Run **independently** of the `mixed_flux=none` test so we learn whether the floor alone
      suffices, the flux alone suffices, or both are needed.
- [ ] `make test` green (the four benchmarks byte-exact; floor unset there).

### Dependencies
Depends on: Phases 1–3.

## Testing Strategy
- Phase 1: config-parse unit (present/absent/negative).
- Phase 2: direct `SolveLSW_TPV205` unit at a tensile σ_n (free-slip vs floored, numbers above);
  TPV205 byte-exact parity (floor off).
- Phase 3: RS strength unit + TPV102/104 byte-exact parity (floor off).
- Phase 4: cluster rerun (bounded V_max) + full `make test`.
- Cross-cutting: a grep test that the floor reaches **all four** strength sites (Table above) — a
  missed site silently leaves a free-slip path.

## Risk Assessment
- **Missed strength site** → a residual free-slip path; mitigate with the all-sites grep test and
  the cluster rerun (if it still blows up, a site was missed).
- **Byte-exactness regression** → the default-disabled sentinel + default-`0.0` `SolveLSW`
  parameter must reproduce current expressions exactly; the parity tests guard this. The
  `SolveLSW_TPV205` signature change touches the oracle kernel — the trailing default is the
  safety.
- **Symptom vs root (carry from §6b):** the floor bounds the friction-mediated feedback but the
  σ_n oscillation source (central-flux / sliver spatial instability) remains; expect a *bounded
  but locally spurious* solution at the slivers. Pair with `mixed_flux=none` + mesh repair.
- **Floor masks legitimate opening** — fine for this pure strike-slip SAFS case (no opening
  expected); scope the option off for normal-fault/bimaterial problems.

## Decisions for /code-review
1. **RESOLVED (2026-05-26, user):** use the **active `μ_eff(δ)`** — `τ_strength = μ_eff(δ)·max(σ_n,
   σ_n_floor)` (C⁰-continuous at the threshold, keeps slip-weakening; "constant in σ_n" below the
   floor but still follows `μ(δ)`). NOT the literal `σ_n_floor·μ_s` hard constant (which would jump
   at the threshold). The plan as written already reflects this — no change needed.
2. **Floor the friction-solver σ_n argument** (`:201/:548`, `:226`-RS) too, or only the explicit
   `strength`? (Affects `Theta`/bracket; v1 = strength only.)
3. **Scalar (this plan) vs per-DOF floor** (depth-dependent cohesion à la TPV31)? v1 scalar.
4. **Also cap σ_n itself / the radiated normal traction**, not just the strength? (Bounding V may
   already bound `[[v_n]]`→σ_n; confirm empirically before adding.)
5. Config key location: `[friction] sigma_n_strength_floor_pa` (top-level, applies to both laws)
   vs per-law block. Plan = top-level.
