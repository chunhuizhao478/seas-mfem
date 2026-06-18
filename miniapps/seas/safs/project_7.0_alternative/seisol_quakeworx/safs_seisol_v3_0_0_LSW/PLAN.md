---
title: "PLAN — safs_seisol_v3_0_0_LSW (CSM depth-dependent stress + stress-tied LSW)"
date: "2026-06-15 — rev 1 (DESIGN)"
---

- Date: 2026-06-15 (rev 1, DESIGN — no input files generated yet)
- Author target: SeisSol single-event dynamic rupture, FL=16 linear slip weakening
- Template: `../safs_seisol_v2_0_0_LSW/` (parameters.par, safs_fault.yaml,
  safs_initial_stress.yaml, safs_material_cvm.yaml/.nc, safs_mesh.puml.h5)
- New stress source: `../toolbox/on_fault_stress_projection_csm/`
  `csm_yhsm2013_stress_on_safs_mesh_*` — the C1 k=2.39, MUSCAL density,
  hydrostatic P_p, linear-sampled depth-dependent on-fault stress field.
- Constraint from the user: the rupture must propagate from the SE end
  (s ~ -18 km) all the way to the NW end (s ~ +266 km).
- Status: DESIGN FOR REVIEW. Building the input files is gated on the
  decisions in Section 8. None of this is verifiable locally — full
  propagation can only be confirmed on Frontera/Expanse (project policy:
  no local full-mesh runs).

---

## 0. Goal and scope

Reproduce the V2 LSW SeisSol case structure, but replace V2's single
CONSTANT regional stress tensor with the depth-dependent, CSM-orientation
stress field we just built, and redesign the linear-slip-weakening (LSW)
friction so a single nucleated rupture runs the full SE -> NW length of
the SAFS multi-strand fault.

SCOPE RULE (carried from the stress-projection plan): every parameter
value must trace to SAF-system / southern-California observations or to
setting-independent mechanics (Mohr-Coulomb friction, energy balance).
No values imported from other earthquakes; Parkfield/SAFOD excluded.

---

## 1. What is adopted from V2 unchanged vs changed

UNCHANGED from `safs_seisol_v2_0_0_LSW`:
- Mesh: `safs_mesh.puml.h5` (same fault, 60,658 DR facets; UTM 11N;
  boundary codes 1=free surface, 3=DR fault, 5=absorbing).
- Material: `safs_material_cvm.yaml` + `safs_material_cvm.nc` (MUSCAL CVM
  via easi !ASAGI). V3 reuses these verbatim.
- parameters.par numerics: FL=16, numflux=godunov, Plasticity=0, CFL=0.5,
  ClusteredLTS=2, FixTimeStep=0.1, Format=10 (no volume output),
  refPointMethod=1 with XRef/YRef/ZRef=(0,-1,0), nucleation ramp s_0=0,
  t_0=1.0, fault &Elementwise output mask, energy + free-surface output.
- Deep barrier concept: lock z in (-20000, -15000) m via mu_s -> 1e6.
- Nucleation concept: time-ramped Gaussian over-stress Tnuc_s at the
  hypocenter (606971, 3707270, -4965.62).

CHANGED in V3:
1. Initial fault stress: from V2's single ConstantMap tensor to the
   SPATIALLY VARYING CSM C1 k=2.39 field (Section 3).
2. Friction: from V2's UNIFORM (mu_s=0.43, mu_d=0.12, d_c=1.8) to
   STRESS-TIED spatially varying mu_s and mu_d (Section 4) — this is
   forced by the new field's mu_app spread; uniform friction provably
   cannot meet the SE->NW constraint (Section 4.1).

---

## 2. The CSM C1 k=2.39 stress field (what V3 adopts)

Source artifact (this repo, regenerated 2026-06-15):
`toolbox/on_fault_stress_projection_csm/csm_yhsm2013_stress_on_safs_mesh_fault_stress.vtu`
+ `_summary.json`. Built with:
`--magnitude-mode depth-R --closure ratio --k-ratio 2.39
 --density-source muscal --pp-model hydrostatic --axes csm --sample linear`.

Per-facet effective tractions (compression-positive, SEAS convention),
fault-wide stats:

```default
sigma_n_eff:  min 0.20   median 111.7   max 415.4   MPa  (grows with depth)
tau_mag:      min 0.06   median  29.2   max 160.2   MPa
mu_app=tau/sn: min 0.014 median  0.275  max  0.449
  percentiles: p75 0.340  p90 0.376  p95 0.397  p99 0.423
  fraction mu_app > 0.40: 4.3 %   > 0.35: 21.9 %
hypocenter facet (s=0): sigma_n_eff 70.1, tau 20.8, mu_app 0.296, z -5054 m
```

ALONG-STRIKE mu_app profile (seismogenic band z = -12..-3 km, strike
azimuth 314 deg, s measured from the hypocenter, NW positive; SE end at
s ~ -18 km, NW end at s ~ +266 km):

```default
 s (km)    mu_app(med)  sigma_n_eff(med MPa)   note
 -20..0      0.246          95.8
   0..20     0.352          81.0               hypocenter (good nucleation)
  20..40     0.262         102.0
  40..60     0.142         120.4               *** PRIMARY GATE (p10 0.103)
  60..80     0.246         125.0
  80..100    0.376         140.9
 100..120    0.268         118.2
 120..140    0.240         123.4
 140..160    0.309         109.2
 160..180    0.352         135.0
 180..200    0.231         184.9
 200..220    0.308         150.0
 220..240    0.181         169.4               *** SECOND GATE (p10 0.144)
 240..260    0.333         174.5
```

The two low-mu_app corridors (s=40-60 km, mu_app ~0.14; s=220-240 km,
mu_app ~0.18) are the propagation gates. The s=40-60 km gate is the
mechanical analog of the San Gorgonio restraining knot: the CSM stress
is poorly oriented for slip on the local fault strike there, so the
resolved shear-to-normal ratio collapses. This gate is DEEPER (lower
mu_app) than V2's ~0.25 corridor where candidates B and C arrested.

---

## 3. Stress injection: CSM field -> SeisSol via easi !ASAGI

V2 fed SeisSol a single global Cartesian tensor (ConstantMap) and let
SeisSol project it per facet (BaseDRInitializer rotateStressToFaultCS).
The CSM field is NOT a single tensor; it is a spatially varying tensor
field sigma_eff(x,y,z). SeisSol's supported mechanism for that is a
gridded easi !ASAGI field (MFEM_to_SeisSol_mapping.md Section 2, line:
"a gridded field via easi !ASAGI (HDF5)"), identical in form to the
material `safs_material_cvm.nc` workflow.

KEY PROPERTY (why this is clean for stress): the CSM effective tensor is
a genuine VOLUMETRIC field — its orientation is the CSM principal-axis
field interpolated in (x,y) (depth-invariant), its shape ratio R is
CSM(x,y), and its magnitude is Sv_eff(z) (MUSCAL lithostat - hydrostatic
P_p) combined with the closure. So sigma_eff(x,y,z) can be evaluated on
any 3D grid and SeisSol projects it per facet exactly as it did the V2
constant tensor.

### 3.1 New tool: CSM stress -> ASAGI nc
Add `toolbox/on_fault_stress_projection_csm/csm_stress_to_asagi.py`
(sibling of `project_csm_stress_to_vtu.py`, reusing its tensor builder):

```default
inputs : the CSM csv (orientations), --k-ratio 2.39, MUSCAL nc (density),
         pp-model hydrostatic   (same knobs as the projection deliverable)
grid   : rectilinear UTM 11N covering the fault bounding box with margin
         x: ~360..625 km, y: ~3690..3840 km, z: -17000..0 m
         spacing: dz = 250 m (match material nc); dx=dy = 500 m
         (the fault facets are ~430 m; 500 m resolves the along-strike
          mu_app gate structure; finer only needed near junctions, see
          Section 4.4)
compute: at every grid node evaluate sigma_eff(x,y,z) (e,n,u = x,y,z),
         then NEGATE to SeisSol compression-NEGATIVE convention.
         P_p is already removed (effective tensor) -> baked-in, exact under
         projection (isotropic part projects to sigma_n only).
output : COARDS NetCDF, compound variable stress{s_xx,s_yy,s_zz,s_xy,
         s_yz,s_xz} laid out data(z,y,x), Pa  (mirror convert_cvm_to_asagi)
selfcheck: trilinear sample of the written nc at the 60,658 fault facet
         centroids must reproduce the projection VTU's per-facet
         sigma_n_eff / tau within tolerance (the projection VTU is ground
         truth). FAIL the converter if not.
```

### 3.2 safs_initial_stress.yaml (V3)
Replace V2's ConstantMap with:
```default
[s_xx,s_yy,s_zz,s_xy,s_yz,s_xz]: !ASAGI
  file: safs_stress_csm.nc
  parameters: [s_xx,s_yy,s_zz,s_xy,s_yz,s_xz]
  var: data
  interpolation: linear
```
(with a !ConstantMap fallback to V2's tensor for graceful degradation,
following the material yaml pattern.)

---

## 4. Friction redesign (the core of V3)

### 4.1 Uniform LSW friction is INFEASIBLE for this field (proof)
To avoid spontaneous failure at t=0, the static friction must exceed the
initial ratio everywhere: mu_s > max(mu_app) = 0.449, so mu_s ~ 0.46.
With a uniform mu_s, the strength excess at the primary gate is

```default
SE_gate = (mu_s - mu_app)*sigma_n = (0.46 - 0.142)*120 = 38.3 MPa
```

The largest stress drop physically available there (even at mu_d = 0) is

```default
drop_max = (mu_app - 0)*sigma_n = 0.142*120 = 17.0 MPa
```

so the S-ratio S = SE/drop >= 38.3/17.0 = 2.25 (and >= 3.6 at the
mu_app ~ 0.10 patches), far above the ~1.5 unbounded-growth threshold.
No uniform (mu_s, mu_d) passes the gate; it would require a negative mu_d.
The user's mu_d ~ 0.10 was tested explicitly: at the gate it yields a
5 MPa drop against a 38 MPa strength excess (S = 7.5), and ZERO drop on
the mu_app ~ 0.10 patches. mu_d is not the binding knob — mu_s is. This
is the same arrest mechanism that stopped V2 candidates B and C, now made
worse by the deeper gate.

### 4.2 Stress-tied friction (the design)
Tie mu_s and mu_d to the local apparent friction so the strength excess
and stress drop (hence S) are CONTROLLED along the entire fault:

```default
mu_s(facet) = mu_app(facet) + delta_s
mu_d(facet) = max( mu_app(facet) - delta_d , mu_d_floor )
=>  strength excess SE = delta_s * sigma_n      (uniform in mu)
    stress drop     dt = delta_d * sigma_n      (where unfloored)
    S = delta_s / delta_d                        (uniform, low)
```

RECOMMENDED values (Section 8 lists these as the decisions to confirm):

```default
delta_s   = 0.05      (t=0 margin; SE = 5% of sigma_n; mu_s = mu_app+0.05)
delta_d   = 0.15      (dynamic drop = 15% of sigma_n where unfloored)
mu_d_floor= 0.00      (full weakening on the lowest-mu_app patches)
=> S = delta_s/delta_d = 0.33 everywhere  (robust, well below 1)
```

Resulting ranges over the fault:
```default
mu_s = mu_app + 0.05 :  0.06 .. 0.50   (mu_s > mu_app by 0.05 => t=0 stable)
mu_d = max(mu_app-0.15, 0):  0 .. 0.30  (MEDIAN ~0.125 ~ the user's ~0.10)
primary gate (mu_app 0.142, sn 120): mu_s 0.192, mu_d 0  (floored),
   drop = 0.142*120 = 17 MPa, SE = 0.05*120 = 6 MPa, S = 0.35  -> PASSES
hypocenter (mu_app 0.296, sn 70):    mu_s 0.346, mu_d 0.146,
   drop = 0.15*70 = 10.5 MPa, SE = 3.5 MPa, S = 0.33
highest mu_app 0.449:                mu_s 0.499, mu_d 0.299  (no t=0 failure)
```

This is a standard heterogeneous-prestress technique: by referencing the
friction to the prestress, the fault is everywhere uniformly (and weakly)
super-critical once nucleated, so a single nucleation runs the full
length. The low absolute mu_d (median ~0.125, floor 0) is consistent with
the SAF weak-fault / low-heat-flow evidence (Lachenbruch and Sass 1980,
1992; Hardebeck and Hauksson 2001; Townend and Zoback 2004) — Section 9.

### 4.3 Cohesive-zone resolution (d_c)
Linear slip-weakening cohesive-zone length (Day et al. 2005 form):
```default
Lambda = (9*pi/32) * mu_shear * d_c / ((mu_s - mu_d) * sigma_n)
```
With (mu_s - mu_d) = delta_s + delta_d = 0.20 (uniform by construction),
mu_shear ~ 24 GPa, d_c = 2.0 m:
```default
sigma_n = 120 MPa (gate)  : Lambda = 1473 m  (~3 facet widths @ 430 m)
sigma_n = 300 MPa (deep)  : Lambda =  589 m  (~1.4 widths -> MARGINAL deep)
sigma_n =  70 MPa (hypo)  : Lambda = 2526 m  (well resolved)
```
RECOMMEND d_c = 2.0 m (V2 used 1.8 m). The deep facets (sigma_n > ~250
MPa, z < ~-12 km, near the locked barrier) are marginally resolved;
mitigations if a convergence check fails: (a) raise d_c, or (b) tie
d_c to sigma_n so Lambda is constant. Deferred until the first run shows
whether the deep band matters (it is adjacent to the mu_s=1e6 barrier).

### 4.4 Implementation of stress-tied friction (the hard part)
mu_app = |tau| / sigma_n is a FAULT-SURFACE quantity — it depends on the
fault normal, so unlike the stress tensor it is NOT a clean volumetric
function of (x,y,z) and cannot be written as a Lua f(x,y,z). It must be
PRE-COMPUTED per facet (the projection tool already does, knowing each
normal) and delivered to SeisSol's fault easi, which is evaluated at
fault Gauss points (x,y,z). Two mechanisms:

```default
Option A (RECOMMENDED): rasterize per-facet mu_s, mu_d to a 3D ASAGI nc.
  The converter computes mu_app per facet, forms mu_s/mu_d per Section 4.2,
  and writes a grid where each node carries the value of the NEAREST fault
  facet (within a band; benign fill outside). SeisSol's linear ASAGI
  lookup at a fault point returns ~ that facet's value (neighbors vary
  smoothly under --sample linear). safs_fault.yaml: [mu_s] / [mu_d] !ASAGI.
  RISK: at multi-strand junctions (strands < grid-cell apart) nearest-facet
  is ambiguous; mitigate with dx=dy=250 m near the fault or a signed band.

Option B: provide per-facet values keyed by an along-fault parameterization
  if SeisSol's fault input supports it (needs verification in the QuakeWorx
  SeisSol build; not assumed here).
```
The same rasterization can also emit fault-local tractions Ts0/Td0/Pn0 if
we prefer per-facet stress over the volumetric tensor of Section 3 — that
guarantees mu_s - mu_app = delta_s EXACTLY (Section 8 decision 2).

### 4.5 Deep barrier and shallow note
- Deep barrier: keep V2's lock, mu_s = 1e6 for z in (-20000, -15000) m,
  via the same Lua override layered on top of the ASAGI mu_s (or baked
  into the rasterized mu_s field).
- Free surface: NOTE the low-mu_app gate is at DEPTH (s=40-60 km,
  z -12..-3 km), NOT near the surface. The shallow band (z > -3 km) has
  mu_app median 0.277 with sigma_n only ~23 MPa, so the stress-tied rule
  there gives mu_d ~ 0.13 and small absolute drops — not a propagation
  bottleneck. No special free-surface mu_d treatment is required; the
  user's "mu_d ~0.1 near the free surface" is already the fault-median
  outcome of Section 4.2.

---

## 5. Nucleation and friction parameters (hypocenter worked calc)
Worked nucleation/friction sheet for the SE-end hypocenter, updated from
the V2 draft to (a) the C1 k=2.39 field facts and (b) the stress-tied
friction of Section 4. Corrections vs the V2 draft are flagged inline.
(delta_s=0.05, delta_d=0.15, d_c=2.0 m are the Section 4.2 / 8 proposed
values, pending confirmation.)

### 5.1 Facts within the hypocenter (V3, C1 k=2.39 field)
Nucleation location (606971, 3707270, -4965.62); nearest DR facet centroid
(607030, 3707289, -5054 m), 108 m away — in the SE region where mu_app ~
0.30-0.35 (the easiest place to nucleate). Stress from the C1 projection
summary; material from `toolbox/hypocenter_facts/hypocenter_facts.py`:

```default
Effective normal stress  sigma_N = 70.1 MPa   (was 63.8: P_p is now
                                  HYDROSTATIC ~49.6 MPa at 5 km, not 20)
Shear stress magnitude    tau_0   = 20.8 MPa   (was 21.2)
Apparent friction         mu_app  = 0.296      (was 0.332)
Shear-wave speed          Vs      = 3032 m/s   (MUSCAL; SeisSol-linear 3002)
Density                   rho     = 2551 kg/m3 (SeisSol-linear 2541)
Shear modulus             mu      = 23.5 GPa   (SeisSol-linear 22.9)
Poisson ratio             nu      = 0.275      (SeisSol-linear 0.278)
```
Material is UNCHANGED from V2 (same MUSCAL point); only the stress changed,
because V3 adopts the C1 field, not V2's constant regional tensor. Note the
C1 hypocenter traction is more strike-parallel (tau_strike 19.8, tau_dip
6.4 -> rake ~ +18 deg) than V2's (rake -58 deg) — favourable for a
right-lateral strike-slip start.

### 5.2 Linear slip-weakening friction (stress-tied; Section 4)
CORRECTION vs the V2 draft: a UNIFORM mu_s = 0.43 is INFEASIBLE for the C1
field — the fault-wide max mu_app is 0.449 (not 0.42), so the 4.3 % of
facets with mu_app > 0.40 would PRE-SLIP at t=0. Requiring mu_s above the
overall apparent friction therefore forces mu_s >= 0.45 if uniform, which
then cannot rupture the s=40-60 km gate (Section 4.1). V3 ties friction to
the local mu_app, so AT THE HYPOCENTER:

```default
mu_s = mu_app + delta_s     = 0.296 + 0.05 = 0.346   (> mu_app: no pre-slip)
mu_d = max(mu_app - delta_d, 0) = 0.296 - 0.15 = 0.146   (fault median ~0.125)
(mu_s - mu_d) = 0.20         d_c = 2.0 m
```

Process zone (mesh must resolve it with 4-5 elements):
```default
Lb = mu*d_c / ((mu_s-mu_d)*sigma_N) = 23.5e9*2.0 / (0.20*70.1e6) = 3352 m
   mesh edge dx ~ 430-500 m (500 m mesh)  ->  Lb/dx ~ 6.7-7.8 elements  OK
   CORRECTION: the V2 draft's "Lb=2138 > 4*dx=2500" was false (2138<2500);
   here Lb=3352 exceeds 4*dx for any dx in [430,625] m.
   Deeper facets shrink Lb (sigma_N grows): Lb=1958 m at the gate
   (sigma_N=120), 783 m near the deep barrier (sigma_N=300) -> marginal
   there only (Section 4.3).
```

Static-friction sanity (no pre-slip):
```default
mu_s = 0.346 > mu_app = 0.296 at the hypocenter, AND mu_s = mu_app + 0.05
everywhere by construction -> the whole fault is below failure at t=0.
```

Strength excess:
```default
SE = tau_s - tau_0 = (mu_s - mu_app)*sigma_N = 0.05*70.1e6 = 3.5 MPa
   (was 6.2 MPa; smaller by design -> nucleation is EASIER under V3.)
```

Nucleation length (Uenishi-Rice 2003; the patch must cover it):
```default
Lnuc = 1.158*mu*d_c / ((1-nu)*(mu_s-mu_d)*sigma_N)
     = 1.158*23.5e9*2.0 / (0.725*0.20*70.1e6) = 5355 m   (half-length 2677 m)
   CORRECTION: the draft's "Lnuc = mu*d_c/(...)" header was the PROCESS-ZONE
   formula; the nucleation length is the U-R form above (with 1.158/(1-nu)).
```

Nucleation patch (gradual over-stress: time-ramped Gaussian Tnuc_s, T_nuc=1 s):
```default
R = 6000 m  (e-fold radius; patch diameter 2R = 12000 m > Lnuc 5355 m)
dtau_nuc(min) = 1.1*tau_s - tau_0 = 1.1*0.346*70.1 - 20.8 = 5.9 MPa
   choose dtau_nuc = 10 MPa  (peak = tau_0 + 10 = 1.27*tau_s; robust margin)
over-stressed footprint r_os = R*sqrt(ln(dtau_nuc/SE))
   = 6000*sqrt(ln(10/3.5)) = 6143 m  >  Lnuc/2 = 2677 m   (onset ratio 2.3)
   CORRECTION: the draft's r_os line mixed R=3500/5500/6000 and SE=6.24;
   with R=6000 and SE=3.5 MPa -> r_os=6143 m. Even dtau=6 MPa gives
   r_os=4399 m > Lnuc/2, so nucleation is robust; 10 MPa keeps V2's margin.
```

Bottom line: under V3 the SE-end hypocenter nucleates MORE easily than V2
(smaller strength excess, larger over-stressed footprint), so R=6000 m +
dtau=10 MPa is ample. The hard problem is PROPAGATION through the
s=40-60 km gate (Section 4.1), which the stress-tied mu_s addresses — not
nucleation. Re-run `hypocenter_facts.py --case-dir safs_seisol_v3_0_0_LSW`
once the V3 stress/friction nc exist to confirm these facts on the actual
case inputs.

---

## 6. Implementation plan (files, once Section 8 is decided)
Phase 1 - stress converter + nc:
  - create `toolbox/.../csm_stress_to_asagi.py` (Section 3.1), run it ->
    `safs_seisol_v3_0_0_LSW/safs_stress_csm.nc`; converter self-check passes.
Phase 2 - friction fields:
  - extend the converter (or a sibling) to emit `safs_friction_csm.nc`
    (mu_s, mu_d per Section 4.2, rasterized per 4.4); self-check vs VTU.
Phase 3 - case files (copy V2, edit):
  - `parameters.par` (V2 verbatim except EndTime; Section 8 decision 5),
  - `safs_initial_stress.yaml` -> !ASAGI safs_stress_csm.nc (Section 3.2),
  - `safs_fault.yaml` -> [mu_s]/[mu_d] !ASAGI safs_friction_csm.nc + deep
    barrier Lua + d_c/cohesion ConstantMap + Tnuc_s Lua (Section 5),
  - reuse `safs_material_cvm.yaml` + `.nc` + `safs_mesh.puml.h5`,
  - `README_port.md` (run guide + the validation gates of Section 7).

---

## 7. Validation gates (cluster; cannot be checked locally)
- [ ] Converter self-checks pass (stress nc and friction nc reproduce the
      projection VTU per-facet values within tolerance) — LOCAL, must pass
      before submitting.
- [ ] t=0 stability: at t=0 NO facet is at failure (max mu_app < min local
      mu_s by construction = delta_s); confirm from the first step's fault
      output Ts0/Td0/Pn0 vs the projection VTU (should match to interp).
- [ ] Nucleation: slip initiates at the hypocenter patch and the cohesive
      zone is resolved (Lambda >= ~3 facets except the deep barrier band).
- [ ] GATE CROSSING: the rupture front passes s=40-60 km (the San Gorgonio
      gate) and s=220-240 km without arrest (rupture-time RT field
      continuous across both); this is THE acceptance criterion.
- [ ] Full propagation SE (s ~ -18 km) -> NW (s ~ +266 km); moment does not
      plateau before the NW end (the failure signature of V2 B/C).
- [ ] EndTime margin: at >= 3 km/s the ~285 km length needs ~95 s; set
      EndTime with margin (Section 8 decision 5).

---

## 8. Open decisions (confirm before the build)
1. Friction model = stress-tied mu_s/mu_d (Section 4.2). CONFIRM, and the
   values: delta_s = 0.05, delta_d = 0.15, mu_d_floor = 0.00 (S = 0.33;
   median mu_d ~ 0.125). Tighter S (smaller delta_d) = stronger drops /
   lower mu_d; looser delta_s = safer t=0 margin but harder gates.
2. Stress mechanism: (A) volumetric tensor ASAGI + SeisSol projects
   (clean, but SeisSol-recomputed mu_app differs from the tool's by
   interpolation, so the delta_s margin is approximate), vs (B) rasterized
   per-facet Ts0/Td0/Pn0 (exact margin, same machinery as the friction nc).
   RECOMMEND B for an exact, guaranteed-stable t=0 margin.
3. Friction injection mechanism: Option A (rasterized ASAGI) vs verifying a
   native per-facet path in the QuakeWorx SeisSol build (Section 4.4).
4. d_c = 2.0 m uniform vs d_c tied to sigma_n for constant Lambda
   (Section 4.3) — affects only the deep band.
5. EndTime (parameters.par): 150-200 s (V2 used 200 s) given ~285 km.
6. Grid resolution of the stress/friction nc: dx=dy=500 m vs 250 m near the
   fault junctions (Section 4.4 junction risk).

---

## 9. References (SAF-system / mechanics; verification status as of 2026-06-15)
- Lachenbruch, A.H. and Sass, J.H. (1980, 1992). Heat-flow bound on
  SAF-resolved shear (~<=20 MPa average) -> low apparent/dynamic friction.
- Hardebeck, J.L. and Hauksson, E. (2001), JGR 106, 21859-21882. Southern
  California: weak SAF in a low-strength crust; crustal deviatoric stress
  ~10 MPa (Landers stress rotation) — supports low mu_d.
- Hardebeck, J.L. and Michael, A.J. (2004), JGR 109, B11303. Intermediate-
  strength SAF (mu ~ half a strong fault), the basis for k=2.39 / mu~0.45
  in the stress field this case adopts.
- Townend, J. and Zoback, M.D. (2004), GRL 31, L15S11; Jones (1988),
  JGR 93, 8869. SHmax ~65 deg to the southern SAF -> weak fault, low
  resolved friction.
- The LSW cohesive-zone-length estimate Lambda (Section 4.3) is the same
  standard (9*pi/32) form V2 used (Day-type / Palmer-Rice; e.g. Day, S.M.
  et al. 2005, JGR 110, B12307) — setting-independent mechanics; verify
  the exact citation/pages before any manuscript use.
- Uenishi-Rice nucleation half-length (Section 5; Uenishi and Rice 2003,
  JGR 108(B1), 2042) — the same nucleation criterion V2 used;
  setting-independent mechanics; verify pages before manuscript use.
