# Reference papers: SAF principal stress magnitudes & directions for simulation input

Compiled 2026-06-09 via adversarially-verified deep research (22 claims
confirmed 2-1 or better by independent verification against the source
abstracts/full texts; 3 over-generalized claims refuted and excluded).
Purpose: justify the regional stress tensor (SHmax, Shmin, Sv, P_p, psi)
prescribed in the SAFS stress pipeline (`stress/code/`, see
`WORKFLOW_fault_stress_z0embed_triq_shmax120_psi69_az23.md`).

Notation: psi = angle between SHmax and the local SAF strike;
mu_app = |tau| / sigma_n_eff on the fault; compression positive.

---

## 1. The verified anchor papers

### Townend & Zoback (2004), GRL 31, L15S11 — doi:10.1029/2003GL018918
**The most simulation-ready quantitative constraint.** [verified 3-0 x3]
- Southern SAF (~400 km): mean **psi = 68 +/- 7 deg** from 70 focal-mechanism
  inversion locations within 10 km of the trace (63 +/- 10 deg using the 34
  points within 5 km). Incompatible with Byerlee-strength fault (needs ~30 deg).
- Assuming critically stressed crust (crustal mu = 0.6, hydrostatic P_p,
  strike-slip regime, evaluated at 7.5 km):
  - psi = 61 deg -> on-fault tau ~ 50 MPa, sigma_n_eff ~ 150 MPa
  - psi = 75 deg -> tau ~ 30 MPa, sigma_n_eff ~ 170 MPa
  - i.e. **mu_app ~ 0.2-0.3** ("moderately weak"), southern CA.
- Central CA (SF Peninsula, psi ~ 85 deg within ~1 km of the fault):
  tau ~ 10 MPa, sigma_n_eff ~ 180 MPa, **mu_app ~ 0.06** (very weak).
- Implied Sv ~ 190 MPa at 7.5 km -> **Sv gradient ~ 25.3 MPa/km**.
- Caveats: tractions are point estimates at 7.5 km (not measured gradients);
  the paper concedes psi < 68 cannot be ruled out within ~5 km of the fault.

### Hickman & Zoback (2004), GRL — doi:10.1029/2004GL020043  [in hand]
SAFOD pilot hole (Parkfield), the source of this repo's hz1671 parameter set.
- **psi is depth-dependent in the only direct deep measurement: 25 +/- 10 deg
  at ~1.0-1.15 km rotating to 69 +/- 14 deg at ~2.05-2.2 km** (~28 deg
  rotation). A constant psi is itself a modeling choice. [verified 3-0]

### Chery, Zoback & Hickman (2004), GRL 31, L15S13 — doi:10.1029/2004GL019521
2.5-D finite-element model fit to far-field + SAFOD pilot-hole stress.
[verified 3-0 x3]
- Only **fault effective friction mu_eff < 0.1** through the seismogenic
  thickness matches both data sets.
- Best model reproduces only ~15 deg of the observed ~28 deg psi rotation.
- Predicts fault-zone interior mean stress **~2x lithostatic** at low shear —
  caveat: do NOT extend the regional tensor unchanged into the fault core.
- Model-dependent (one FE model), constraint depth ~0.8-2.2 km.

### Yang & Hauksson (2013), GJI 194, 100-117 — doi:10.1093/gji/ggt113
SCEC Community Stress Model contribution (1981-2010 SCSN focal mechanisms).
[verified 3-0 x2]
- **Regional SHmax azimuth ~ N7E averaged over southern California** —
  usable directly; psi then varies along strike as the SAF strike rotates
  (Big Bend). Documents real lateral heterogeneity (NNW SHmax near
  Cajon/Tejon Pass, NNE near San Jacinto/ECSZ): N7E is an average only.
- Reproduces Hardebeck & Hauksson (1999): psi decreases from ~60 deg far
  field to **~40 deg within 10-20 km of the fault** — the contested
  intermediate-angle alternative (damped inversions of Hardebeck & Michael
  2006 reduce the apparent rotation).

### Hardebeck & Michael (2004), JGR 109, B11303 — doi:10.1029/2004JB003239
The intermediate-angle alternative. [verified 2-1 x2 — live dispute]
- Their compiled orientations are NOT consistent with the classic
  weak-fault-in-strong-crust model; observed psi ~ 40-65 deg.
- Alternative (a): **intermediate-strength SAF** — mu = 0.29-0.54,
  fault strength 65-85 MPa at 10 km (130-165 MPa transpressional),
  crustal deviatoric ~100 MPa.
- Alternative (b): **all major faults weak** — strength ~ earthquake stress
  drop (~10 MPa), crustal deviatoric ~ tens of MPa, psi varying 45-90 deg
  through the seismic cycle.

### Zoback & Healy (1992), JGR 97, 5039 — doi:10.1029/91JB02175  (Cajon Pass)
The strong-crust / weak-fault framework. [verified 3-0 x2]
- In-situ magnitudes at 0.9-3.5 km: crust adjacent to the SAF is at
  frictional failure with **Byerlee mu 0.6-1.0 under hydrostatic P_p**.
- SAF itself slides at shear stress comparable to seismic stress drops
  (a few MPa) — consistent with the missing frictional heat-flow anomaly.
- Sets the off-fault stress level + the hydrostatic-P_p assumption; caps
  on-fault initial shear traction in weak-fault models.

### Shamir & Zoback (1992), JGR 97 — doi:10.1029/91JB02959  (Cajon Pass)
- Apparent SHmax = **N57E +/- 19 deg** (32,616 breakout determinations,
  1.75-3.46 km, 4.2 km from the SAF). Against the local ~N60W strike:
  zero resolved RIGHT-lateral shear (actually left-lateral sense) on
  SAF-parallel planes. [verified 3-0 x2]
- Strongly depth-variable (breakout rotations up to ~100 deg); possibly a
  local Cleghorn-fault anomaly (Scholz & Saucier 1993). Calibrates the
  weak-fault geometry; not usable alone as a regional psi.

### McGarr, Zoback & Hanks (1982), JGR 87, 7797 — doi:10.1029/JB087iB09p07797
Off-fault deviatoric magnitude gradient (Mojave, 29 hydrofrac measurements,
80-849 m, 2-34 km from the SAF). [verified 3-0 x2, 2-1]
- Maximum shear stress depth gradient **~7.9 MPa/km** ("typical for
  continents"); extrapolates to ~56 MPa average shear in the upper 14 km.
- If regional deviatoric stress were limited by SAF strength: implied fault
  **mu ~ 0.45** (intermediate). Use as an off-fault deviatoric gradient,
  NOT as on-fault resolved shear in a weak-SAF model.

### Zoback & Roller (1979), Science 206, 445
- Shear stress increases with distance from the SAF (fault-perpendicular
  profile near Palmdale, shallow <0.85 km); on-fault shear inferred to grow
  slowly with depth, reaching only stress-drop magnitudes. [verified 2-1/3-0]
- Original observational basis for a weak-fault on-fault traction gradient.

### Hickman, Zoback & Healy (1988), JGR 93, 15183 (Hi Vista, Mojave)
[verified 3-0 x3]
- At <0.6 km depth, 32 km from the SAF: regime **transitional between
  thrust and strike-slip (SHmax > Shmin ~ Sv)** — pure strike-slip
  Andersonian ordering is wrong at shallow depth.
- SHmax azimuth scatter NNE-NW: unusable as a psi input.
- Magnitudes at 32 km <= those 4 km from the fault: weakens the
  lateral-heterogeneity constraint at depth.

---

## 2. Candidate parameter sets for simulation input

Common scaffolding (all sets): Sv = lithostatic ~25-26 MPa/km; hydrostatic
P_p ~ 9.8-10 MPa/km off-fault; vertical SAF; Andersonian strike-slip at
seismogenic depth (transitional SHmax > Shmin ~ Sv above ~0.6 km);
SHmax/Shmin magnitudes backed out from the critically-stressed-crust
condition (crustal mu = 0.6 on optimally oriented planes).

| Set | psi (deg) | on-fault tau | sigma_n_eff | mu_app | basis |
|---|---|---|---|---|---|
| A. Weak fault, central CA | ~85 (63-85 near-fault) | ~10 MPa @7.5 km | ~180 MPa @7.5 km | ~0.06 | Townend & Zoback 2004 (SF Peninsula); Zoback & Healy 1992; Chery et al. 2004 |
| B. Moderately weak, southern CA | 61-75 (mean 68+/-7) | ~30-50 MPa @7.5 km | ~150-170 MPa @7.5 km | ~0.2-0.3 | Townend & Zoback 2004; SHmax az ~N7E (Yang & Hauksson 2013) |
| C. Intermediate-strength SAF | 40-65 | 65-85 MPa @10 km (SS); 130-165 transpress. | from critical crust, ~100 MPa deviatoric | 0.29-0.54 (~0.45 McGarr 1982) | Hardebeck & Michael 2004; McGarr et al. 1982 |
| D. All-faults-weak | 45-90 (cycle-dep.) | ~10 MPa (~stress drop) | low crustal deviatoric (~tens MPa) | ~0.05-0.1 | Hardebeck & Michael 2004 alt. model |

Depth dependence: in A/B/D on-fault tau grows slowly with depth (order of
stress drops) while off-fault deviatoric stress grows at ~7.9 MPa/km;
near-fault psi is depth-dependent in the only direct deep measurement
(25 -> 69 deg over 1.0-2.2 km, SAFOD pilot hole).

---

## 3. Where this repo's current parameter set sits

Current run (`..._triq_shmax120_psi69_az23`): SHmax=120, Shmin=55, Sv=50,
P_p=20 MPa, psi=69 deg, constant with depth.

- **psi = 69 deg == the Townend & Zoback southern-SAF mean (68 +/- 7)** —
  squarely Set B. The psi-sweep tail-area optimum (72.5 deg, see
  `results/psi_sweep_triq_shmax120_pp20/`) is inside T&Z's 61-75 bracket.
- The magnitudes correspond to an anchor depth of ~2 km on the verified
  gradients: Sv = 50 MPa at 25.3 MPa/km -> z ~ 1.98 km; hydrostatic P_p
  there ~ 19.6 ~ 20 MPa.
- (SHmax - P_p)/(Sv - P_p) = 100/30 = 3.33 -> crustal frictional state
  mu ~ 0.64: a critically stressed crust at the bottom of the Byerlee
  0.6-1.0 range (Zoback & Healy 1992). This is exactly why the psi sweep
  found the field-max mu_apparent pinned at 0.639 for all psi: the
  orientation bound (s1'-s3')/(2 sqrt(s1' s3')) IS the critical-crust mu.
- Shmin = 55 slightly > Sv = 50: the transitional SS/thrust ordering seen
  at shallow Mojave depths (Hickman et al. 1988) — consistent with a
  ~2 km anchor, but note the regime should rotate to SHmax > Sv > Shmin
  at seismogenic depth, which a constant-with-depth tensor cannot do.

## 4. Coverage gaps & unverified leads

The verification stage killed or never reached two requested topics — treat
the following as leads, not verified references:

- **SCEC CSM products / stress-magnitude models**: SCEC CSM portal
  (southern.scec.org/research/csm); GJI 211(1):472 (2017) (likely
  Luttrell & Smith-Konter, crustal differential stress limits from
  topography + focal mechanisms); JGR doi:10.1029/2020JB020817.
- **How dynamic-rupture/SEAS papers parameterize pre-stress**: Geosphere
  18(6):1710 (2022) "The effects of pre-stress assumptions on dynamic
  rupture..."; Science Advances doi:10.1126/sciadv.1500621; SCEC SEAS BP5
  spec (strike.scec.org/cvws/seas/download/SEAS_BP5.pdf).

Other verified-debate context: Lachenbruch & Sass heat flow (Science 238:
1105, 1987 framing); Scholz (2000, Geology) strong-fault dissent and
Scholz & Saucier (1993) Cajon Pass reinterpretation; Lockner et al. (2011)
SAFOD core mu ~ 0.15 saponite gouge; Faulkner et al. (2006) damage-zone
stress rotation alternative; Rice (1992) elevated fault-zone pore pressure
(fault-core P_p is unconstrained by any verified claim here).

## 5. Open questions for this project

1. Near-fault psi: 40-65 deg (Hardebeck school) vs 63-85 deg (Zoback
   school) is unresolved; our sweep shows the critically-stressed-area
   metric is fairly flat over 68-75 deg, so Set B is robust to this within
   the Zoback-school bracket, but Set C (psi ~ 50) would change the picture.
2. Depth dependence: constant-with-depth magnitudes anchor the tensor at
   ~2 km; a lithostatic_sv run (gradients are already supported by
   `project_to_fault_stress.py`) would be needed for seismogenic-depth
   consistency (Sv_grad ~ 25-26 MPa/km, P_p_grad ~ 9.8-10 MPa/km).
3. Fault-zone interior: regional tensor should not be blended unchanged
   into the fault core (Chery et al. 2004 ~2x lithostatic mean stress;
   Rice 1992 overpressure).
