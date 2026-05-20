# Spatial dynamic rupture — post-event debug (job 7738686, 2026-05-20)

First long SAFS dynamic-rupture run (`spatial_dyn_smoke`, `tfinal = 20 s`,
projected H&Z stress, μ_s=0.65 / μ_d=0.50, 1000 m lcfar3000 mesh, 400 ranks).
Prior runs were `tfinal = 0.5 s` and never reached the post-nucleation regime,
so both symptoms below are newly exposed, not regressions.

Two distinct symptoms, **both rooted in the spatial heterogeneity of the
apparent friction** `μ_app(x) = |τ_pre|/σ_n_eff` across the fault combined
with mesh quality — see analysis.

---

## Symptom

### Issue [1] — nucleation does not propagate spontaneously
- The `gradual_overstress` patch nucleates a small event (V_max peaks
  ~0.66 m/s at t≈0.9 s, ParaView shows a red slip-rate patch near the
  hypocenter), but the rupture **does not spread**: the rest of the fault
  stays dark blue (locked).  V_max decays back to ~0.02 m/s by t≈2 s.
- User hypothesis: under-resolved mesh (`L_nuc/h_min = 9.49 < 10`).

### Issue [2] — spurious velocity blowup at the north tip after t≈3 s
- A localized grid-scale speckle at the **shallow north tip** generates
  `slip_rate_strike` that explodes: V_max 0.025 → 9.4e9 → 1.3e34 → … →
  4.2e130 → NaN (then V_max reports 0 because `NaN > x` is false).
- Picked points (ParaView, fault.vtkhdf): id=63 (448252, 3.7984e6, **−100**),
  id=67 (448338, 3.79844e6, **−249.7**) — consistently the shallow north end.
- The TPV-param slots in the picked point (`fault_x2/x3`, `param_a/Dc`) read
  0 — that is **expected** (the SAFS driver leaves those empty; the static
  SAFS params live in `volume.vtkhdf`), not a bug.

---

## Evidence gathered

### A. Friction regime varies across the fault (from config + run log)
`--print-derived` and the projected-stress config give, at the hypocenter
(σ_n_eff ≈ 71.5 MPa, τ_pre ≈ 32.5 MPa) vs the north tip (μ_app,max = 0.6410,
documented in the config header at (477917, 3780180, −181.6)):

| location | μ_app = τ_pre/σ_n_eff | static ratio τ_pre/(μ_s·σ_n) | dynamic stress drop τ_pre − μ_d·σ_n |
|---|---|---|---|
| hypocenter | 0.455 | 0.70 | **32.5 − 35.75 = −3.25 MPa  (NEGATIVE)** |
| north tip  | 0.641 | **0.986** | positive, and τ_pre/(μ_d·σ_n) ≈ **1.28** |

With a single uniform μ_s=0.65 / μ_d=0.50, `μ_app(x)` spans **0.45 → 0.64**,
which **straddles μ_d = 0.50**.  That single fact drives both symptoms.

### B. Mesh quality at the north tip (local 1000 m mesh, fault = phys 101)
Parsed `safs_fault_box_nwcut_1000m_lcfar3000.msh` (10 603 fault triangles)
and scored triangle quality (min interior angle, aspect, 2·r_in/r_circ):

```
ALL fault triangles:   min_angle median 58.5°, p10 52°, min 6.4°
                       aspect    median 1.03,  p90 1.17, max 7.92
                       radius_ratio median 0.999

Closest triangles to (448338, 3.79844e6, −249.7):
  d=324 m  min_edge 177 m  aspect 5.78  min_angle  9.89°  rr 0.308  z=−150
  d=407 m  min_edge 177 m  aspect 6.34  min_angle  8.67°  rr 0.259  z=−150
  d=652 m  min_edge 243 m  aspect 4.09  min_angle 14.08°  rr 0.430  z=−223

5 worst slivers MESH-WIDE (min_angle):
   6.40°  aspect 7.92  (429883, 3807432, −141)
   7.35°  aspect 7.31  (429271, 3807727, −141)
   8.67°  aspect 6.34  (447972, 3798588, −150)
   9.43°  aspect 3.19  (365053, 3838902, −1903)
   9.89°  aspect 5.78  (448576, 3798244, −150)
```
**The worst sliver triangles in the entire fault mesh are clustered at the
shallow north fault-surface trace (z = −141…−150 m), exactly where the
blowup originates.**  Median fault triangles are near-equilateral (58°), so
this is a localized meshing defect, not a global one.

The project triangle-q metric (`4√3·A/Σedge²`, 1 = equilateral) at the tip
slivers is **0.26–0.29** (mesh-wide worst 0.19).  NOTE: the canonical
`meshing/code/check_mesh_quality.py` gate is `min eta/q > 0.1`, so these tip
slivers **pass the existing gate** — the gate is too lenient to catch them
for dynamic rupture.  Reproduce:
  - spatial localisation (pure stdlib, no deps):
    `meshing/code/locate_fault_slivers.py MESH.msh --target 448338 3798440 -249.7`
  - canonical η/q gates + histograms (needs `meshio`; run under
    `conda activate pythonenv`):  `meshing/code/check_mesh_quality.py MESH.msh`

### C. Code findings (read-only)
- Dynamic flux uses `strength = |σ_n_total| · f_V` (fault_face_flux.cpp:209)
  → tensile σ_n over-resists rather than opening; not a blowup source.
- `friction_solver.cpp:80` guards `σ_n ≤ 0 → V = τ/η`.
- `FaultGeometry` has degenerate-normal / t1-degenerate fallbacks
  (fault_geometry.hpp:1206, R-103/R-304).  The run logged
  `zero-normal fallback: 0`, but a **sliver** triangle does NOT trip the
  zero-normal guard — it yields a *computable but ill-conditioned* normal /
  face quadrature, which is the more likely failure mode here.
- `min_sigma_n_pa = 1 MPa` floor is applied (spatial_friction.cpp:1036), but
  depth_model="constant" means σ_n is set by fault **orientation**, not
  depth — so the shallow north tip is high-μ_app because of geometry, not a
  σ_n-floor artifact.

---

## Root-cause analysis

### Issue [1]: the nucleation sits in a NEGATIVE-stress-drop region
For LSW, spontaneous propagation requires a non-negative dynamic stress drop
`τ_pre − μ_d·σ_n ≥ 0` along the path (positive energy release).  At the
hypocenter the stress drop is **−3.25 MPa** (τ_pre 32.5 < μ_d·σ_n 35.75).
So once the `gradual_overstress` forcing reaches full ramp and the rupture
tries to leave the forced patch, it enters background fault that **absorbs**
rather than releases energy → it arrests.  The observed event is essentially
just the forced patch slipping.

**This is a physics/calibration problem, not (primarily) a resolution
problem.**  A finer mesh cannot make a negative-stress-drop fault propagate.
The user's `L_nuc/h_min = 9.49` under-resolution is a *real secondary* issue
(it degrades the cohesive zone and would impair propagation even with a
positive stress drop), but it is not the dominant cause here.

The only dynamically supercritical region (μ_app > μ_d = 0.5) is the **north
end** (μ_app up to 0.64) — which is also the over-critical, sliver-meshed
region that blows up (Issue 2).  So with this calibration there is **no clean
"nucleate-and-propagate-a-controlled-event" regime**.

### Issue [2]: sliver elements at the over-critical north tip → grid-scale instability
Three factors coincide at (≈448000, 3.798e6, −150):
1. **Mesh:** the worst sliver triangles in the fault (min_angle 6.4–9.9°,
   aspect ~6–8) — ill-conditioned face normals and DG face integrals →
   numerical noise injected into the local traction / slip-rate.
2. **Stress:** the most over-critical point on the fault — static ratio
   0.986 (≈0.4 MPa from yield) and dynamic ratio 1.28 (supercritical: once
   weakened it cannot arrest).
3. **Timing:** after the main event, afterslip + radiated/reflected waves
   perturb this near-critical point past yield (~t≈3 s).

The slivers **seed** grid-scale noise; the over-criticality **amplifies** it
(no arrest mechanism); the explicit ADER scheme grows it ~30%/step → NaN.
The base scheme is otherwise stable (2000 clean steps + a clean event), so
this is a localized seeded instability, not a global CFL failure.

### Q: will reducing the time step help?  (Likely NOT — and here is why)
`WaveOperator::ComputeMaxDt` sets the CFL length scale to the tet **inscribed
diameter** `6·Vol/Σface_area` (wave_operator.inl:61-63), *not* the edge
length.  So the global dt **already accounts for sliver geometry** — a sliver
tet has a small inscribed diameter and drives dt down.  Consistency check:
the run's `dt_cfl = 0.0015 s` ⇒ effective h ≈ dt·cp/cfl ≈ 18 m, far below the
231 m edge `h_min`, i.e. dt is already set by the worst sliver's inradius.

Therefore the north-tip slivers are **already CFL-covered**; the blowup is
not a CFL violation but (a) the *geometric ill-conditioning* of the sliver
fault-face normal / quadrature (independent of dt) and (b) the over-critical
physical runaway (also dt-independent).  Halving dt (`SAFS_CFL=0.25`) is
worth running as the discriminator, but the expectation is the blowup
re-appears at the **same physical time (~3 s)**.  If so, dt is exonerated and
the fixes are the mesh (slivers) and the calibration (criticality), not the
time step.  (Also: a global dt cut to chase a few tip slivers is very
expensive — 20 s already ≈ 13k steps at dt=0.0015 s.)

---

## Recommended actions (in priority order)

1. **Fix the mesh at the fault-surface trace (addresses Issue 2 seed).**
   Re-mesh so the fault∩free-surface intersection at the north end is not
   slivered — e.g., refine the surface trace, add a conforming buffer, or
   improve the nwcut booleanization there.  Target min_angle ≳ 25–30° and
   aspect ≲ 3 on all fault triangles.  Re-run `check_mesh_quality.py` to
   confirm the slivers are gone before any production run.

2. **Re-examine the friction/prestress calibration (addresses both).**
   `μ_app(x) ∈ [0.45, 0.64]` straddling μ_d = 0.50 is the core problem.
   Options (modeling decision — needs the user):
   - Lower μ_d so the nucleation region has a positive stress drop
     (μ_d < 0.45 at the hypocenter) — but this makes the north *more*
     supercritical (worse blowup) unless the north mesh/stress is also fixed.
   - Re-site the nucleation into a region where μ_app > μ_d (so it can
     propagate) that is *not* slivered/over-critical.
   - Reduce the regional-stress shear/normal ratio so μ_app,max drops well
     below μ_s and the fault is uniformly subcritical-but-nucleatable.

3. **PML / shorter tfinal (mitigates the t≈3 s perturbation), secondary.**
   No PML today; reflections from the 1st-order absorbing boundary reach the
   fault around t≈3 s and help tip the near-critical north point.  `--pml`
   or `tfinal ≤ ~3 s` removes that trigger but does NOT fix the underlying
   sliver + over-criticality.

## Discriminating experiments still open (from the A/B knobs already wired)
- `SAFS_CFL=0.25` (halve dt): if the north blowup time is unchanged → it is
  physical/seeded (mesh+stress), not time stepping.  Expected: unchanged.
- μ_s/μ_d shift (keep Δμ=0.15 so L_nuc/resolution is unchanged): tests the
  criticality contribution independent of resolution.
- Mesh refinement at the fault trace: the decisive test for Issue 2 — if the
  blowup disappears with healthy north-tip elements (same stress), the
  slivers were the seed.

## Status / what is confirmed vs. inferred
- CONFIRMED (data): north-tip slivers (mesh parse); negative hypocenter
  stress drop (config numbers); the |σ_n| flux handling and fallback logic
  (code read); the blowup location and growth (run log + ParaView).
- INFERRED (needs a re-run to confirm): that healthy north-tip elements
  remove the blowup; that the negative stress drop (not resolution) is the
  dominant cause of non-propagation.  Both are testable with the knobs above.
