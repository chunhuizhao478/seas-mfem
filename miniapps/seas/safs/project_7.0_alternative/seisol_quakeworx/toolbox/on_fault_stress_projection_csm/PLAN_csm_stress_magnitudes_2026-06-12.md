---
title: "PLAN — Depth-dependent stress magnitudes for the CSM (YHSM-2013) fault projection"
date: "2026-06-15 — rev 7"
---

- Date: 2026-06-15 (rev 7)
- Status: PLAN (no implementation yet; current tool implements constant
  prescribed principals 80/40/35 MPa).  Four of the open decisions are
  now FIXED by user directive 2026-06-15 (density = MUSCAL, pore
  pressure = hydrostatic, closure = Luttrell differential stress,
  axes/taper unchanged) — see Section 8.
- Tool: `project_csm_stress_to_vtu.py` (this folder)
- Data (orientations): `raw_data/yang_and_hauksson/CSM_data_1781290479811.csv`
  (YHSM-2013, orientation only)
- Data (differential stress, NEW rev 5):
  `raw_data/luttrell_differential_stress/CSM_data_1781537253009.csv`
  (Luttrell-2017 CSM export; D = S1-S3 column, see Section 3 closure C2-s)
- Mesh: `safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5`
  (UTM 11N, fault z from -21 m to -16,528 m)
- Revision history:
  - rev 7 (2026-06-15) adds the `--sample {nearest,linear}` orientation
    sampler (Section 5). Nearest-neighbour sampling of the ~2 km CSM
    grid made the mu / tau field blocky (60,658 facets -> ~719 distinct
    CSM cells); `--sample linear` barycentrically interpolates the CSM
    tensor components per facet for a smooth orientation field (NN
    fallback outside the data hull). Default stays `nearest`
    (byte-identical to prior runs); the mu cap is preserved under either.
  - rev 6 (2026-06-15) replaces the closure-C2 positivity cap. The
    rev-5 implementation capped sig3 to a fixed fraction of Sv_eff,
    which pinned the stress ratio sig1/sig3 at ~30 on shallow facets
    (super-critical) and produced unphysical apparent friction mu -> 1.0
    at the free surface. The cap is now a FRICTIONAL (Mohr-Coulomb /
    Byerlee) ceiling: D_applied = min(D, D_max) with
    D_max = (k_max-1)*Sv_eff/((1-R)*k_max+R), bounding sig1/sig3 <= k_max
    (CLI `--cap-k-max`, default 5.83 = mu 1.0) and tapering the
    differential to 0 as Sv_eff -> 0 at the surface. Adds references 23
    (Sibson 1974), 24 (Brace & Kohlstedt 1980), 25 (Townend & Zoback
    2000); see Sections 3 (C2-s cap), 4, 7. The seismogenic-depth field
    and the hypocenter are unchanged (those facets were never capped).
  - rev 5 (2026-06-15) records four user decisions and the arrival of
    the Luttrell spatial differential-stress data: (1) DENSITY = the
    MUSCAL lateral-mean rho(z) profile is the fixed default; the
    Vp->rho-recipe confirmation is downgraded from a blocker to a
    manuscript footnote (Section 2.1, decisions 5+7). (2) PORE PRESSURE
    = hydrostatic is the fixed default (Section 2.2, decision 2).
    (3) CLOSURE = the differential stress D = S1-S3 from Luttrell &
    Smith-Konter (2017) is now the recommended default, used SPATIALLY
    (closure C2-s) from the CSM-portal export that has been obtained and
    placed in `raw_data/luttrell_differential_stress/` — this resolves
    decisions 1 and 6 and flips the rev-3 recommendation away from the
    frictional ratio C1 (Sections 0, 3, 4, 6). (4) The C1-vs-C2
    comparison is expanded into a dedicated subsection (Section 3.x) and
    the worked example is rerun across depths with the on-fault Luttrell
    field (Section 4). The frictional ratio C1 is retained as the
    `--closure ratio` alternative.
  - rev 4 corrects the density provenance: the project's material file
    derives from the raw extracts in
    `velocity/raw/multiscale_statewise_cvm/` — the SCEC CVM Explorer
    model "muscal" (Multi-scale Statewide CALifornia model, Yeh &
    Ben-Zion, CANVAS base) — NOT CVM-H, and rev 3's claim that the
    densities are Nafe-Drake/Brocher conversions is downgraded to
    "conventional expectation, to be confirmed against the muscal
    publication" (Section 2.1).
  - rev 3 addresses four review comments: (1) verified the CVM density
    variable and clarified bulk-vs-dry density for the lithostatic
    integral (Section 2.1); (2) replaced the bare Hubbert-Rubey
    justification with SAF-specific fluid-pressure measurements and
    models (Section 2.2); (3) removed all Parkfield/SAFOD-derived
    values — creeping central SAF, not representative of the locked
    southern SAF (Sections 3, 6, 7); (4) added a path to using the
    Luttrell & Smith-Konter constraint as a SPATIAL field
    (Section 3, closure C2-s).
  - rev 2 removed every reference to, and every parameter value taken
    from, the Kaikoura (New Zealand) dynamic rupture literature
    (Ulrich et al. 2019); those are NOT San Andreas constraints.

---

## 0. Goal, approach, and scope rule

The CSM model YHSM-2013 (Yang & Hauksson 2013) provides stress
ORIENTATION only — its tensors are normalized deviatoric ("This model
provides orientation only. Magnitude information is not meaningful.",
csv header). The current tool therefore prescribes three constant
effective principal magnitudes. This plan upgrades the magnitudes to a
depth-dependent, observation-constrained model in three steps:

1. S2 = effective vertical stress Sv_eff(z) = total overburden
   (lithostatic, BULK density from MUSCAL) minus hydrostatic pore
   pressure — Step 1.
2. A literature closure linking the two horizontal principals S1, S3.
   The DEFAULT (rev 5) is the differential stress D = S1 - S3 taken as
   a SPATIAL field D(x,y) from Luttrell & Smith-Konter (2017) — closure
   C2-s. The frictional effective ratio k = S1/S3 (closure C1) is kept
   as an alternative — Step 2.
3. Combine Sv_eff(z), the closure, and the PER-POINT stress shape ratio
   R from the CSM csv (column 13) to solve for all three magnitudes
   S1(x,z) >= S2(z) >= S3(x,z) in closed form — Step 3.

Then rebuild the regional tensor from the CSM angles (principal-axis
eigenvectors / SHmax azimuth, csv columns 10, 21–35) and project onto
the fault with the existing machinery — Steps 4–5.

**SCOPE RULE (standing constraint for this work).** This is a San
Andreas Fault model. Every parameter VALUE must trace to
(a) observations in the SAF system / southern California, or
(b) setting-independent mechanics (effective-stress principle,
Mohr-Coulomb friction, Anderson faulting theory). Two exclusions:

- Values calibrated for OTHER earthquakes or regions (Kaikoura,
  Landers-event specifics, subduction settings, ...) must not be
  imported, however similar the methodology looks.
- Parkfield / SAFOD data (pilot-hole stress, SAFOD core gouge
  friction) are EXCLUDED from value-carrying roles: Parkfield sits in
  the creeping central SAF, mechanically unlike the locked southern
  SAF this mesh models (user directive 2026-06-12).

Methodology-only citations are allowed but must be flagged as such.

---

## 1. Conventions — READ FIRST (sign/ratio traps)

All formulas below use SEAS compression-positive effective principal
stresses, written

```default
sig1 >= sig2 >= sig3 > 0     (compression positive, EFFECTIVE, MPa)
```

The csv is TENSION-positive ("Tension is taken as positive"), with
S1 = most tensional. Mapping between the two:

```default
sig1 = -S3_csv   (csv most COMPRESSIONAL principal -> our largest)
sig2 = -S2_csv
sig3 = -S1_csv   (csv most TENSIONAL principal -> our smallest)
```

Shape-ratio dictionary (verified numerically on all 22,131 csv rows,
agreement to 5e-4 = csv rounding):

```default
csv col 13:  R_GF  = (S2-S3)/(S1-S3)   tension-positive
                   = (sig1-sig2)/(sig1-sig3)   compression-positive

csv col 12:  phi   = (S1-S2)/(S1-S3)   tension-positive  (Angelier 1979)
                   = (sig2-sig3)/(sig1-sig3)   compression-positive
                   (this compression-positive form is often written
                    "nu" or "Phi" in the dynamic-rupture literature)

identity:    phi + R_GF = 1   (holds exactly in the csv)
```

So: when this plan says "R from the csv" it means column 13. Do not
mix the two columns.

Useful interpolation forms (both equivalent):

```default
sig2 = (1 - R)*sig1 + R*sig3            with R   = csv col 13
sig2 = phi*sig1 + (1 - phi)*sig3        with phi = csv col 12
```

Frame: (east, north, up) identified with mesh (x, y, z); UTM Zone 11N
metres; lon/lat -> UTM via pyproj EPSG:4326 -> EPSG:32611, always_xy
(velocity-pipeline convention). Grid convergence (< ~1 deg) neglected.

---

## 2. Step 1 — S2 = effective vertical stress Sv_eff(z)

```default
Sv_total(z) = integral_0^z rho_bulk(z') * g * dz'    (total overburden)
P_p(z)      = pore pressure model (Section 2.2)
Sv_eff(z)   = Sv_total(z) - P_p(z)
```

### 2.1 Density: what is in the MUSCAL material file, and bulk vs dry

**Do we have density?** Yes — verified 2026-06-12.
`safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc` stores a compound
variable `data{rho, mu, lambda}` on a 264 x 195 x 181 UTM grid
(dz = 250 m, z = 0 .. -45 km, units attribute "rho kg/m^3"). The
`rho` field is carried through unmodified from the raw extracts in
`velocity/raw/multiscale_statewise_cvm/` (`Lon,Lat,Vp(m/s),Vs(m/s),
Density(kg/m^3)` columns; see `velocity/code/raw_readers.py` and
`toolbox/generate_velocity_nc_from_raw`).

**Which velocity model is this?** The raw slice headers say
`# CVM(abbr): muscal` — the SCEC CVM Explorer export of MUSCAL, the
Multi-scale Statewide CALifornia velocity model (Yeh, Ben-Zion &
Olsen, ESSOAr preprint 2026; UCVM >= 25.7 component
`SCECcode/muscal`, authors file: Te-Yang Yeh and Yehuda Ben-Zion,
USC). MUSCAL builds on the CANVAS adjoint-waveform-tomography base
model of California and Nevada (Doody et al. 2023) and merges
regional/local datasets with a near-surface low-velocity taper. Note
this is NOT CVM-H; the `velocity/README.md` header line predates the
switch (three raw versions exist — `cvmh/`, `cvm_s4.26.m01/`,
`multiscale_statewise_cvm/` — and the production sidecar recorded in
the nc attribute `source_sidecar` is the muscal one).

Computed lateral-mean profile (this mesh region):

```default
z =      0 m   mean rho = 1678 kg/m^3   (basin sediments at surface)
z =  -1000 m   mean rho = 2450
z =  -2000 m   mean rho = 2525
z =  -5000 m   mean rho = 2669
z = -10000 m   mean rho = 2773
z = -15000 m   mean rho = 2759
global range 1245 .. 3341 kg/m^3
```

**Is it rock density or rock+fluid combined?** It is in-situ BULK
density (rock matrix + pore fluid), with one caveat on provenance.
The `rho` array ships inside the authors' own model file
(`model_MUSCAL_CANVAS_dll0.01_vardz_float32_cmpd.nc`; the UCVM
packaging only re-encodes it — verified in `SCECcode/muscal`
`data/rewrite2bin.py`). Seismic tomography does not invert for
density, so MUSCAL's rho is necessarily ASSIGNED from the velocities;
the universal practice is an empirical Vp-density relation
(Nafe-Drake curve as formulated by Brocher 2005), whose calibration
data are saturated, in-situ rocks — i.e. bulk density. The exact
relation MUSCAL used is not stated in the UCVM packaging. The shallow
values (lateral mean 1678 kg/m^3 at the surface, minimum 1245) are
only consistent with saturated-sediment bulk densities, supporting
this reading.

DECISION (user directive 2026-06-15, resolves decisions 5+7): USE the
MUSCAL density. The rho field shipped in `safs_material_cvm.nc` is
taken AS-IS — its lateral-mean rho(z) profile (table below) is the
fixed default for the lithostatic integral. The Vp->rho-recipe
confirmation that rev 3/4 listed as a pre-implementation ACTION is
DOWNGRADED to a non-blocking manuscript footnote: cite the bulk-density
reading as the conventional expectation (Brocher 2005), to be replaced
by the model's own documentation when the Yeh, Ben-Zion & Olsen
preprint or software@scec.org confirms the exact relation. It does not
gate implementation.

**And that is the CORRECT density for this integral.** The total
vertical stress at depth is the weight of EVERYTHING above —
mineral grains and the fluid filling their pores:

```default
Sv_total(z) = integral rho_bulk * g * dz      (bulk = correct)
```

Using a dry-skeleton density would UNDERESTIMATE the overburden. The
"dry rock minus fluid" intuition from the original task statement is
realized by the effective-stress subtraction Sv_eff = Sv_total - P_p
(Hubbert & Rubey 1959), NOT by integrating a dry density. So:
bulk density in the integral, pore pressure subtracted after — no
double counting.

Magnitude check (lateral-mean MUSCAL profile vs constant 2670 kg/m^3):

```default
z = -5.0  km: Sv_cvm = 123.6 MPa   const-2670 = 131.0   (-5.7 %)
z = -10.0 km: Sv_cvm = 258.3 MPa   const-2670 = 261.9   (-1.4 %)
z = -16.5 km: Sv_cvm = 434.4 MPa   const-2670 = 432.2   (+0.5 %)
```

The difference is largest in the upper 5 km (low-density basin
sediments — Coachella/Salton region), exactly where the fault mesh is
shallow; this is why the MUSCAL profile is the DEFAULT over a constant.
Implement as `--density-source {muscal,constant}` (default `muscal`)
with `--material-nc PATH` (the production file is
`safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc`; the name retains the
"cvm" tag but the contents are MUSCAL, see above). `constant` (with
`--rho-const-kgm3`, default 2670) is kept only for the magnitude check
and reproducibility. Use the laterally-AVERAGED rho(z) profile, not
per-column integration: per-column Sv would introduce lateral Sv jumps
that are not in mechanical equilibrium and would leak into mu maps.

### 2.2 Pore pressure: SAF-specific observations and models

What Step 1 needs is P_p(z) in the REGIONAL field (country rock) —
the fault-core fluid state is a separate, local matter. The available
SAF-system constraints (searched 2026-06-12; there is NO community
pore-pressure model for southern California analogous to CVM/CFM/CSM):

- DIRECT measurement, southern California, 4 km from the SAF:
  Cajon Pass scientific borehole. In situ fluid-pressure build-up
  tests at ~1.8–2.1 km depth give pressures only ~5 % above
  hydrostatic (Coyle & Zoback 1988). Independently, the measured
  stress magnitudes over 0.9–3.5 km are consistent with frictional
  equilibrium at mu = 0.6–1.0 UNDER HYDROSTATIC pore pressure
  (Zoback & Healy 1992).
- SAF-specific fluid-pressure DISTRIBUTION models: Rice (1992)
  (idealized) and Fulton & Saffer (2009) (numerical, fed by the
  mantle-helium fluid source of Kennedy et al. 1997) both conclude
  that overpressure develops INSIDE the low-permeability fault zone
  while the surrounding crust remains near-hydrostatic. Fulton &
  Saffer quantify the fault-zone EXCESS pressure: ~1 MPa at 2.7 km,
  ~17 MPa at 6 km depth, confined to the fault zone and mostly below
  5 km.

DECISION (user directive 2026-06-15, resolves decision 2): USE
hydrostatic pore pressure for the regional field. This is the
defensible model and is the fixed default,

```default
P_p(z) = rho_w * g * z          (rho_w = 1000 kg/m^3)
lambda(z) = P_p/Sv_total = 0.37 - 0.40 for the MUSCAL profile above
            (0.397 at 5 km, 0.380 at 10 km, 0.373 at 16.5 km)
```

CLI: `--pp-model hydrostatic` (default). A FAULT-LOCAL effective
normal stress reduction following the Fulton & Saffer (2009) excess-
pressure depth profile is a possible later refinement — it would be
applied on the fault surface after projection, not to the regional
tensor; out of scope for this tool, noted for the friction setup.

(Hubbert & Rubey 1959 is cited only for the effective-stress
principle and the lambda notation — methodology, not values.)

### 2.3 Optional deep taper of the deviatoric part

Below the seismogenic zone the crust creeps and cannot sustain large
deviatoric stress, so the deviatoric part should taper out with depth
while the isotropic effective part continues. SAF-specific depth
scale: the southern California seismogenic thickness is 15.0 (+1.2/
-1.1) km on average, locally < 10 km in the Salton Trough and > 20 km
elsewhere (Nazareth & Hauksson 2004) — directly relevant since the
SAFS mesh spans the Coachella/San Gorgonio corridor and reaches
z = -16.5 km.

```default
sigma_eff(x,z) = Omega(z) * sigma_dev_eff(x,z) + Sv_iso_eff(z) * I
Omega(z) = 1 above z_seis, smooth -> 0 below     (default OFF)
--taper-zseis-km 15      (Nazareth & Hauksson 2004 regional average)
```

Keep OFF by default until first dynamic-rupture trial; the fault mesh
bottom (16.5 km) is near the regional seismogenic floor anyway.

### 2.4 Regime validity check (S2 = vertical assumption)

Assigning the vertical stress to S2 assumes a strike-slip (Andersonian)
regime, Aphi in [1, 2] (Anderson faulting classes; Aphi of Simpson
1997). Measured ON THIS FAULT MESH by nearest-neighbour sampling of
the csv (computed 2026-06-12):

```default
Aphi < 1 (normal-transtensional)   16.8 % of facets
1 <= Aphi <= 2 (strike-slip)       76.2 %
Aphi > 2 (reverse-transpressive)    7.0 %    (San Gorgonio corridor)

V2 (intermediate axis) plunge:  median 73.6 deg
  plunge >= 60 deg: 67.8 %   plunge >= 45 deg: 76.2 %
```

So S2 ~ vertical is a good approximation on ~3/4 of the fault and
questionable on the rest (consistent with the transpressive San
Gorgonio knot in the regional studies, Yang & Hauksson 2013;
Hardebeck & Hauksson 2001). Mitigation: output the sampled V2 plunge
and Aphi as cell fields (`csm_V2_plunge_deg_cell`, `csm_Aphi_cell`) so
the violation zones are visible in ParaView, and report the
sigma_zz_eff vs Sv_eff(z) discrepancy statistics in the summary
(Section 6).

---

## 3. Step 2 — Literature closure between S1 and S3 (SAF region)

One more scalar relation is needed to close the system. Two usable
forms, both to be implemented. ALL supporting values below are from
the SAF system / southern California; Parkfield/SAFOD data are
excluded per the scope rule.

### Closure C1 — effective stress ratio k = sig1/sig3 (ALTERNATIVE, rev 5)

Frictional-equilibrium limit on optimally oriented faults
(Mohr-Coulomb; friction range mu = 0.6–1.0 of Byerlee 1978):

```default
k = sig1/sig3 = ( sqrt(mu^2 + 1) + mu )^2
mu = 0.6  ->  k = 3.12
mu = 1.0  ->  k = 5.83
```

SAF-region evidence:

- Cajon Pass borehole, 4.2 km from the SAF, southern California
  (Zoback & Healy 1992): measured shear/normal stress ratios on
  favorably oriented planes at 0.9–3.5 km depth match mu = 0.6–1.0
  with hydrostatic pore pressure. This is the closest direct
  magnitude measurement to the modeled fault section and the primary
  justification for k ~ 3 (mu = 0.6, lower bound of the measured
  range) in the REGIONAL field.
- Townend & Zoback (2004): stress fields adjacent to the SAF in
  central and southern California are consistent with a STRONG crust
  while the SAF itself slips at high SHmax angle (up to ~85 deg) —
  "strong crust, weak fault". The CSM ORIENTATIONS supply the
  misorientation; k supplies the strong-crust magnitude.
- Strong-fault counterpoint: Scholz (2000) argues the heat-flow
  evidence is not decisive and the SAF may sustain Byerlee-level
  stress. Either way k ~ 3 is the defensible regional value; the
  weak-vs-strong debate concerns the resolved stress ON the fault,
  which here comes from the CSM orientations.

### Closure C2 — differential stress D = sig1 - sig3 (DEFAULT, rev 5)

Luttrell & Smith-Konter (2017), from southern California topography +
focal mechanisms: differential stress at seismogenic depth must
EXCEED ~20 MPa in most of southern California and ~62 MPa near the
most rugged SAF topography (which includes San Gorgonio, on this
mesh). As of rev 5 their model has been obtained as a gridded field
(closure C2-s below), so D is used SPATIALLY rather than as a single
dialed scalar. Use D as a direct dial:

```default
D options: spatial D(x,y) from the Luttrell field (DEFAULT, C2-s);
           or a constant scalar D (e.g. 30-60 MPa, --diff-MPa);
           or D proportional to Sv_eff(z) (self-similar with depth;
           equivalent to a fixed k).
```

A note on depth dependence. The Luttrell field is published at a
single reference depth (5 km below sea level) and the authors state
explicitly that "these model values do not vary with depth" (csv
header). So in closure C2 the differential stress D(x) is held
CONSTANT with depth while the intermediate principal sig2 = Sv_eff(z)
grows with depth — i.e. the deviatoric part is BOUNDED and the
overburden carries the depth dependence. This is the defining contrast
with closure C1 (where D grows linearly with depth); see Section 3.x.

### Closure C2-s — SPATIAL Luttrell & Smith-Konter field (OBTAINED, rev 5)

Review question (rev 3): "can we get the spatial distribution of
this?" Answer (rev 5): YES — the data is now in hand. The Luttrell &
Smith-Konter (2017) model is on the SAME SCEC Community Stress Model
portal that this project's YHSM-2013 csv came from, exported in the
identical `CSM_data_*.csv` 35-column format, and has been placed at:

```default
raw_data/luttrell_differential_stress/CSM_data_1781537253009.csv
  header: CSM model "Luttrell-2017", contributed by K. Luttrell and
          B. Smith-Konter; "provides both meaningful orientation and
          meaningful magnitude"; tension-positive MPa, WGS84 lon/lat.
  19,100 valid rows on a ~0.02 deg (~2 km) lon/lat grid
          (lon -119.149 .. -114.918, lat 32.625 .. 35.252).
  single depth = 5 km below sea level; "values do not vary with depth".
  column 16 (1-based) = "differential stress (S1-S3) in MPa".  <-- USE
  (column 15 "isotropic pressure" is EMPTY for every row: this is a
   deviatoric / differential-stress model, no absolute mean stress.)
```

The grid FULLY contains the SAFS fault footprint (fault lon
-118.474 .. -115.686, lat 33.365 .. 34.684), so every facet has a
nearby sample. Sampling the field over the fault footprint
(2026-06-15, this folder):

```default
on-fault D = S1-S3 (MPa):  min ~0.1   p25 ~12   median ~19   mean ~20
                           p75 ~26    p90 ~34    max ~60
  hypocenter facet (Coachella/Salton end):  D ~ 12 MPa
  San Gorgonio knot (rugged SAF topography): up to ~60 MPa
```

These reproduce the paper's headline bounds (>= ~20 MPa regionally,
~62 MPa at the most rugged SAF topography) directly on the mesh.

Use of D(x) (= column 16) with a SECOND nearest-neighbour sampler,
identical in form to the CSM-orientation sampler (depth ignored — the
field is depth-invariant by construction):

```default
read:    parse the same way as the YHSM csv (lon, lat, col 16 -> D);
         drop rows with NaN D.
sample:  build a cKDTree on the Luttrell (x,y) in UTM 11N; query at
         each fault-facet centroid (x,y); D_facet = D at nearest node.
apply:   closure C2 (Section 4) with this per-facet D, depth-constant:
             sig1 = Sv_eff(z) + R(x)*D(x)
             sig3 = Sv_eff(z) - (1-R(x))*D(x)
cap:     FRICTIONAL ceiling (rev 6). The brittle crust cannot sustain a
         differential stress above the Mohr-Coulomb/Byerlee frictional
         strength, so the stress ratio sig1/sig3 is bounded by
         k_max = (sqrt(mu^2+1)+mu)^2; expressed as a differential cap
         this is the closure-C1 value at k_max:
             D_max(z) = (k_max - 1)*Sv_eff(z) / ((1-R)*k_max + R)
             D_applied = min(D(x), D_max(z))
         Because D_max is PROPORTIONAL to Sv_eff it -> 0 at the free
         surface, so the near-surface differential vanishes instead of
         producing a super-critical ratio. Count and report capped
         facets. (Byerlee 1978; Sibson 1974; Brace & Kohlstedt 1980;
         observed in the SAF region by Zoback & Healy 1992 and Townend &
         Zoback 2000.)
CLI:     --closure diff  --diff-csv PATH  --cap-k-max 5.83  (PATH
         defaults to the Luttrell file above; --diff-MPa is the scalar
         fallback when no csv is given; --cap-k-max default 5.83 = mu 1.0
         Byerlee ceiling, use 3.12 for mu 0.6).
```

Validation against the published bound becomes a self-consistency
check: D_applied(x) == D_Luttrell(x) on every UN-capped facet, and the
capped-facet count quantifies where the frictional ceiling bit (the
shallow + high-D San Gorgonio facets).

### 3.x — C1 vs C2: how the two closures differ (and why C2 is default)

Both closures take sig2 = Sv_eff(z) as the (sub-)vertical intermediate
principal and use the SAME per-point shape ratio R to place sig2
between sig1 and sig3. They differ ONLY in what sets the SPREAD of
sig1 and sig3 around sig2 — and that single difference propagates into
opposite depth behaviour and opposite stress levels.

1) What sets the magnitude scale.

```default
C1 (ratio):  sig1 = k * sig3,  k = ((mu^2+1)^0.5 + mu)^2 fixed scalar.
             The scale is tied to FRICTION (Mohr-Coulomb equilibrium on
             an optimally oriented fault).  No magnitude data needed.
C2 (diff):   D = sig1 - sig3 prescribed directly (here = Luttrell(x)).
             The scale is tied to an OBSERVATIONAL differential-stress
             bound (topography + focal mechanisms).  Needs a D field.
```

2) Depth behaviour — the key divergence. Because sig3 (C1) is
proportional to Sv_eff(z), the C1 differential stress GROWS with depth;
the C2 differential stress is HELD CONSTANT (Luttrell is depth-
invariant). Writing both as a function of Sv_eff:

```default
C1:  D_C1(z) = sig1 - sig3 = (k-1) * sig3
             = (k-1) / ((1-R)*k + R) * Sv_eff(z)     ~ proportional to z
     max shear (sig1-sig3)/2 grows without bound with depth.
C2:  D_C2    = Luttrell(x)  = constant in z
     max shear = D/2 is the SAME at every depth (a few tens of MPa).
```

3) Numbers on this mesh (MUSCAL rho, hydrostatic P_p, on-fault median
R = 0.44; C1 with k = 3.12; C2 with the on-fault median Luttrell
D = 19 MPa):

```default
depth   Sv_eff   C1 max-shear (k=3.12)   C2 max-shear (D=19)
 5 km    74.6 MPa      36.1 MPa                 9.5 MPa
10 km   160.2 MPa      77.7 MPa                 9.5 MPa
16.5 km 272.5 MPa     132.0 MPa                 9.5 MPa
```

C1 climbs to >100 MPa of shear at the base of the seismogenic zone
("strong crust", Byerlee-bounded). C2 stays flat at D/2 — ~9.5 MPa for
a median facet, up to ~30 MPa for the San Gorgonio high-D (D~60)
facets — i.e. at or below the ~20 MPa-class heat-flow / weak-fault
bound (Lachenbruch & Sass 1980, 1992) almost everywhere.

4) Consequences for the rupture model.

```default
- Resolved shear tau and mu_apparent scale with the deviatoric part,
  so C1 produces MUCH higher tau / mu_apparent than C2 (the constant-
  magnitude run, ~C1 at 4-5 km, already sits at mu_app median 0.34).
  C2/Luttrell brings mu_apparent toward the 0.1-0.3 weak-fault target.
- C1 needs no spatial data but cannot represent lateral contrasts;
  C2-s carries the San Gorgonio differential-stress high directly from
  the data (~60 MPa there vs ~12 MPa at the Coachella/Salton hypocenter).
- C2's positivity limit D < Sv_eff/(1-R) bites at shallow depth where
  D is large (capping needed); C1 is automatically ordered for k>=1.
```

5) Why C2/Luttrell is the rev-5 default. The modelled fault is the
LOCKED southern SAF, where the heat-flow and stress-rotation
observations (Sections 3 consistency check, 6) point to a weak fault
in a moderately stressed crust — a BOUNDED deviatoric stress of a few
tens of MPa, exactly what the Luttrell field supplies, and spatially
resolved. C1 (frictional, strong-crust) is retained as
`--closure ratio` for sensitivity studies and the strong-crust
counterpoint (Scholz 2000).

### Consistency check against SAF weak-fault observations

Whatever closure is chosen, the RESOLVED stress on the SAF must be
checked after projection (Section 6) against:

- heat-flow bound: fault-parallel shear stress averaged over the
  seismogenic layer <= ~20 MPa (Lachenbruch & Sass 1980, SAF-wide;
  Lachenbruch & Sass 1992 specifically at Cajon Pass; synthesis
  Zoback et al. 1987);
- low apparent friction on the SAF (~0.1–0.3) implied by the high
  SHmax-to-fault angle in the southern California inversions
  (Hardebeck & Hauksson 2001; Townend & Zoback 2004).

If the projected tau / mu_apparent grossly exceed these, reduce k (or
D) toward the weak end — with the reduction documented against these
same references.

Recommended default (rev 5, user directive): C2-s — the SPATIAL
Luttrell differential-stress field. CLI:

```default
DEFAULT:  --closure diff   --diff-csv <LUTTRELL_CSV>   (spatial, C2-s)
fallback: --closure diff   --diff-MPa 40               (scalar, no csv)
alt:      --closure ratio  --k-ratio 3.12              (frictional C1)
  <LUTTRELL_CSV> = raw_data/luttrell_differential_stress/
                   CSM_data_1781537253009.csv
```

---

## 4. Step 3 — Closed-form magnitudes from (Sv_eff, closure, R)

Given sig2 = Sv_eff(z) and the csv R = R_GF (column 13) at the sampled
CSM point:

### With C1 (ratio k):

```default
sig2 = (1-R)*sig1 + R*sig3   and   sig1 = k*sig3
=>  sig3(x,z) = Sv_eff(z) / ( (1-R(x))*k + R(x) )
    sig1(x,z) = k * sig3(x,z)
```

Ordering is automatic for k >= 1 and R in [0,1]:
sig1 >= sig2  <=>  R*(k-1) >= 0;  sig2 >= sig3  <=>  (1-R)*(k-1) >= 0.

### With C2 (differential D) — DEFAULT, rev 5:

```default
sig1(x,z) = Sv_eff(z) + R(x)*D(x)
sig3(x,z) = Sv_eff(z) - (1-R(x))*D(x)
D(x) = Luttrell(x) (closure C2-s, depth-constant), or scalar --diff-MPa
frictional cap:  D_applied = min(D, D_max),
                 D_max(z) = (k_max-1)*Sv_eff(z) / ((1-R)*k_max + R)
                 => sig1/sig3 <= k_max  (Byerlee ceiling, default 5.83)
```

The naive positivity bound D < Sv_eff/(1-R) is necessary but NOT
sufficient: it keeps sig3 > 0 yet still admits a near-zero sig3, i.e. a
super-critical sig1/sig3 ratio that produces unphysical apparent
friction (mu -> 1.0) at the free surface. The FRICTIONAL cap above
(rev 6) is the physical one: the brittle crust cannot exceed
Mohr-Coulomb/Byerlee frictional strength, so sig1/sig3 <= k_max, and
since D_max is proportional to Sv_eff the differential vanishes as the
overburden -> 0 at the surface (Byerlee 1978; Sibson 1974; Brace &
Kohlstedt 1980; Zoback & Healy 1992; Townend & Zoback 2000). This bites
at shallow depth and where the Luttrell D is large (San Gorgonio high);
count and report capped facets.

### Worked example (on-fault median R = 0.44, MUSCAL profile, hydrostatic P_p):

```default
                      z = 5 km    z = 10 km   z = 16.5 km
Sv_total              123.6 MPa   258.3 MPa   434.4 MPa   (MUSCAL, Sec 2.1)
P_p (hydrostatic)      49.1 MPa    98.1 MPa   161.9 MPa
Sv_eff = sig2          74.6 MPa   160.2 MPa   272.5 MPa

C1, k = 3.12 (frictional, strong crust):
  sig3                 34.1 MPa    73.2 MPa   124.6 MPa
  sig1                106.3 MPa   228.5 MPa   388.7 MPa
  max shear            36.1 MPa    77.7 MPa   132.0 MPa  <- grows with z

C2-s, D = Luttrell (default); median on-fault D = 19 MPa, depth-const:
  sig1                 82.9 MPa   168.6 MPa   280.9 MPa
  sig3                 63.9 MPa   149.6 MPa   261.9 MPa
  max shear             9.5 MPa     9.5 MPa     9.5 MPa  <- flat (= D/2)
  (San Gorgonio facet, D ~ 60 MPa: max shear ~30 MPa, capped if shallow)
```

The depth divergence is the whole point: C1 max shear climbs 36 ->
132 MPa over the seismogenic column, while C2/Luttrell stays flat at
D/2. For comparison the current constant-magnitude run uses 80/40/35
MPa at ALL depths — roughly the C1 numbers at z ~ 4–5 km. The depth-
dependent model removes exactly this inconsistency, and the C2/Luttrell
default keeps the resolved shear within the SAF weak-fault bounds
(Section 3 consistency check).

---

## 5. Step 4 — Assemble the tensor from the csv angles

Per CSM point we have the principal axes; two assembly options:

### Option A (recommended default; matches current tool + user intent)

Use the full eigenvectors of the csv tensor (columns 21–29, or
equivalently our own eigh of See..Suu — already verified identical,
|dot| median 1.0000):

```default
sigma_eff(x,z) = sig1 * u_c u_c^T + sig2 * u_i u_i^T + sig3 * u_t u_t^T
u_c = most-compressional axis (csv V3)
u_i = intermediate axis        (csv V2)
u_t = most-tensional axis      (csv V1)
```

Honest caveat: sig2 = Sv_eff(z) is assigned to u_i, which is only
SUB-vertical (median plunge 73.6 deg on the fault); the built tensor's
sigma_zz_eff will deviate from Sv_eff(z) where the axes tilt. Report
the deviation stats; flag low-plunge zones via cell fields. The
benefit: the measured transpressive/transtensional style (San
Gorgonio vs Salton) is retained.

### Option B (strict Andersonian)

Anderson faulting theory (Anderson 1951): one principal axis exactly
vertical. Use only the SHmax azimuth (csv column 10):

```default
e_H = (sin(az), cos(az), 0)    az = SHmax deg E of N
e_h = (cos(az), -sin(az), 0)   (horizontal, perpendicular)
e_v = (0, 0, 1)
sigma_eff = sig1 * e_H e_H^T + sig3 * e_h e_h^T + sig2 * e_v e_v^T
```

Cleaner mechanics (sigma_zz_eff == Sv_eff exactly), but discards the
measured axis plunges — in the 24% non-strike-slip zones the local
style (e.g. the San Gorgonio thrust component) is lost. Implement
both: `--axes {csm,andersonian}`.

Either way the per-point tensor then goes through the EXISTING
pipeline: CSM-orientation sampling at facet centroids (now also using
the facet centroid DEPTH for Sv_eff), Tandem basis (s = up x n,
d = s x n), normal harmonisation to the SW half-space, traction
resolution, mu_apparent = tau_mag / sigma_n_eff.

Orientation sampling — nearest vs linear (rev 7). The CSM (YHSM-2013)
grid is ~2 km (0.02 deg); the fault mesh is ~430 m triangles, so
nearest-neighbour sampling makes the orientation (and shape ratio R)
PIECEWISE-CONSTANT in the ~2 km CSM cells: on this mesh 60,658 facets
map to ~719 distinct CSM points (median ~78 facets per cell), which
renders as blocky mu / tau patches. `--sample linear` instead does a
barycentric (Delaunay) interpolation of the SIX tension-positive tensor
components onto each facet centroid and re-eigendecomposes per facet to
get a SMOOTH axis + R field (SHmax / eig-gap / Aphi / V2-plunge
diagnostics are recomputed from the interpolated tensor to match);
facets outside the CSM convex hull fall back to nearest neighbour (zero
on this fault — it is fully inside the CSM coverage).

```default
--sample {nearest,linear}   (default nearest, reproducible/blocky;
                             linear = smooth orientation field)
```

Caveats. (1) Smoothing the magnitude scale is unaffected: Sv_eff is
already depth-interpolated and R only sets the intermediate axis;
`--sample` changes only the ORIENTATION continuity. (2) The mu cap is
preserved under either sampling: for closure C1 the ratio k bounds
mu <= (k-1)/(2 sqrt k) for ANY orientation, so interpolation only
redistributes mu within [0, that bound]. (3) Linear interpolation
smooths ACROSS the San Gorgonio style transition (transtension <->
transpression); it interpolates the tensor components (not the discrete
Aphi regime), and the `csm_Aphi_cell` / `csm_eig_gap_cell` outputs flag
where a smoothed facet straddles a style boundary.

---

## 6. Step 5 — Implementation and validation plan

### Implementation (single file change, project_csm_stress_to_vtu.py)

1. Add magnitude model: `--magnitude-mode {constant, depth-R}`.
   `constant` keeps today's behavior (sigma1/2/3 flags) for
   reproducibility of the existing artifact.
2. `depth-R` mode flags, with rev-5 defaults shown:
   `--closure diff` (default; `ratio` for C1),
   `--diff-csv PATH` (default = the Luttrell file in
   `raw_data/luttrell_differential_stress/`; closure C2-s),
   `--diff-MPa 40` (scalar fallback when no csv),
   `--k-ratio 3.12` (only used by `--closure ratio`),
   `--pp-model hydrostatic` (default),
   `--density-source muscal` (default; `constant` + `--rho-const-kgm3`
   for the check) with `--material-nc PATH`
   (default `safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc`),
   `--axes csm` (default; `andersonian` alt),
   `--sample nearest` (default; `linear` = interpolated/smooth
   orientation, see Section 5),
   optional `--taper-zseis-km` (off by default).
3. New per-cell outputs: sig1/sig2/sig3_eff_MPa_cell, csm_R_cell,
   csm_Aphi_cell, csm_V2_plunge_deg_cell, Sv_total/Sv_eff/P_p_MPa_cell,
   depth_m_cell, and (closure diff) D_diff_MPa_cell (the sampled
   differential, column 16 for C2-s), D_applied_MPa_cell (after the
   frictional cap), and a per-facet capped flag; record the
   capped-facet count in the summary.
4. Summary JSON: record every parameter + closure formulas used.

### Validation checklist (all must pass before handing fields to easi)

- [ ] Shape recovery: eigenvalues of every built tensor reproduce the
      input R to ~1e-3 (closed-form, so failure = code bug).
- [ ] Ordering + positivity: sig1 >= sig2 >= sig3 > 0 on 100% of
      facets (count capped facets in C2 mode).
- [ ] SHmax of built tensor vs csv column 10 (median < 1 deg, as in
      the constant-magnitude run).
- [ ] Option A only: sigma_zz_eff vs Sv_eff(z) deviation stats;
      correlate with V2 plunge.
- [ ] SAF weak-fault sanity: distribution of resolved tau on the
      fault vs the ~20 MPa-class heat-flow bound (Lachenbruch & Sass
      1980, 1992); mu_apparent median in the ~0.1–0.3 range implied
      by the high SHmax-to-fault angle (Hardebeck & Hauksson 2001;
      Townend & Zoback 2004). If violated, revisit k / D per
      Section 3.
- [ ] Luttrell self-consistency (closure C2-s, data now in hand):
      D_applied(x) == D_luttrell(x) on every un-capped facet; report
      the capped-facet count and the depths at which capping occurred
      (expected: shallow + San Gorgonio high-D facets). The Luttrell
      nn-sampling distance must stay within the grid spacing (~2 km),
      like the YHSM orientation sampler.
- [ ] Hypocenter facet report (same probe point as previous runs:
      606971, 3707270, -4965.62).
- [ ] Diff against the constant-magnitude artifact
      (`csm_yhsm2013_stress_on_safs_mesh_*`): same orientation fields,
      magnitudes now depth-dependent.

---

## 7. References (exact citations; verification status as of 2026-06-12)

Group 1 — SAF system / southern California observations (these carry
ALL parameter values). Verified online (publisher page) unless noted:

1. Yang, W. & Hauksson, E. (2013). The tectonic crustal stress field
   and style of faulting along the Pacific North America Plate
   boundary in Southern California. Geophysical Journal International
   194(1), 100–117. doi:10.1093/gji/ggt113.
   [Source model YHSM-2013: orientations, phi/R/Aphi, SHmax.]
2. Zoback, M.D. & Healy, J.H. (1992). In situ stress measurements to
   3.5 km depth in the Cajon Pass scientific research borehole:
   implications for the mechanics of crustal faulting. Journal of
   Geophysical Research 97(B4), 5039–5057. doi:10.1029/91JB02175.
   [Borehole 4.2 km from the southern SAF: mu = 0.6–1.0 frictional
   equilibrium with HYDROSTATIC pore pressure -> k and lambda values.]
3. Coyle, B.J. & Zoback, M.D. (1988). In situ permeability and fluid
   pressure measurements at ~2 km depth in the Cajon Pass research
   well. Geophysical Research Letters 15(9), 1029–1032.
   doi:10.1029/GL015i009p01029.
   [DIRECT fluid-pressure measurement near the southern SAF: ~5 %
   above hydrostatic — the key datum behind --pp-model hydrostatic.]
4. Fulton, P.M. & Saffer, D.M. (2009). Potential role of
   mantle-derived fluids in weakening the San Andreas Fault. Journal
   of Geophysical Research 114, B07408. doi:10.1029/2008JB006087.
   [SAF-specific pore-pressure DISTRIBUTION model: overpressure
   localized in the fault zone (~1 MPa @ 2.7 km, ~17 MPa @ 6 km
   excess), country rock near-hydrostatic; supports hydrostatic
   regional lambda + optional fault-local correction.]
5. Townend, J. & Zoback, M.D. (2004). Regional tectonic stress near
   the San Andreas fault in central and southern California.
   Geophysical Research Letters 31, L15S11. doi:10.1029/2003GL018918.
   [Strong crust / weak SAF; SHmax at high angle to the fault.]
6. Hardebeck, J.L. & Hauksson, E. (2001). Crustal stress field in
   southern California and its implications for fault mechanics.
   Journal of Geophysical Research 106(B10), 21859–21882.
   doi:10.1029/2001JB000292.
   [Stress rotations near the SAF; low apparent fault friction;
   mu_apparent validation target.]
7. Luttrell, K. & Smith-Konter, B. (2017). Limits on crustal
   differential stress in southern California from topography and
   earthquake focal mechanisms. Geophysical Journal International
   211(1), 472–482. doi:10.1093/gji/ggx301. (SCEC contribution #6064.)
   [Closure C2/C2-s, the rev-5 DEFAULT: D >= ~20 MPa regionally, ~62
   MPa near rugged SAF topography. The gridded model has been OBTAINED
   from the SCEC CSM portal as `CSM_data_1781537253009.csv` (CSM model
   "Luttrell-2017", placed in `raw_data/luttrell_differential_stress/`):
   35-column tension-positive export, differential stress in column 16,
   single depth 5 km, depth-invariant, ~0.02 deg grid covering the
   whole SAFS fault footprint. Used as the spatial D(x,y) field.]
8. Lachenbruch, A.H. & Sass, J.H. (1980). Heat flow and energetics of
   the San Andreas fault zone. Journal of Geophysical Research
   85(B11), 6185–6222. doi:10.1029/JB085iB11p06185.
   [Heat-flow bound: SAF-resolved shear <= ~20 MPa; validation.]
9. Lachenbruch, A.H. & Sass, J.H. (1992). Heat flow from Cajon Pass,
   fault strength, and tectonic implications. Journal of Geophysical
   Research 97(B4), 4995–5030. doi:10.1029/91JB01506.
   [Same bound specifically at Cajon Pass, i.e. the southern SAF;
   pages re-verify before manuscript use.]
10. Zoback, M.D. et al. (1987). New evidence on the state of stress of
    the San Andreas fault system. Science 238, 1105–1111.
    [Weak-SAF / strong-crust synthesis.]
11. Scholz, C.H. (2000). Evidence for a strong San Andreas fault.
    Geology 28(2), 163–166.
    doi:10.1130/0091-7613(2000)28<163:EFASSA>2.0.CO;2.
    [Strong-fault counterpoint; bounds the debate the validation
    checklist navigates.]
12. Nazareth, J.J. & Hauksson, E. (2004). The seismogenic thickness of
    the southern California crust. Bulletin of the Seismological
    Society of America 94(3), 940–960.
    [z_seis = 15.0 +1.2/-1.1 km average, < 10 km Salton Trough; deep
    taper depth scale.]
13. Kennedy, B.M. et al. (1997). Mantle fluids in the San Andreas
    fault system, California. Science 278, 1278–1281.
    [Mantle-helium evidence feeding the Fulton & Saffer model;
    standard citation, re-verify pages before manuscript use.]

Group 2 — setting-independent mechanics and data relations (no
event-specific values):

14. Hubbert, M.K. & Rubey, W.W. (1959). Role of fluid pressure in
    mechanics of overthrust faulting, I. GSA Bulletin 70, 115–166.
    [Effective-stress principle + lambda notation only; verified.]
15. Anderson, E.M. (1951). The Dynamics of Faulting and Dyke Formation
    with Applications to Britain, 2nd ed., Oliver & Boyd. [Andersonian
    regimes: one principal axis vertical; classic monograph.]
16. Byerlee, J. (1978). Friction of rocks. Pure and Applied Geophysics
    116, 615–626. [Laboratory mu = 0.6–1.0 range quoted by the
    borehole studies; standard citation, re-verify pages before
    manuscript use.]
17. Brocher, T.M. (2005). Empirical relations between elastic
    wavespeeds and density in the Earth's crust. Bulletin of the
    Seismological Society of America 95(6), 2081–2092.
    doi:10.1785/0120050077. [The conventional empirical Vp->rho
    relation for crustal models, calibrated on saturated in-situ
    rocks (-> bulk density). Whether MUSCAL used exactly this
    relation is TO CONFIRM, Section 2.1; citation itself verified.]

Group 3 — definitions used by the csv columns (cited by the csv header
itself; standard citations, re-verify volume/pages before manuscript
use):

18. Gephart, J.W. & Forsyth, D.W. (1984). An improved method for
    determining the regional stress tensor using earthquake focal
    mechanism data... Journal of Geophysical Research 89(B11),
    9305–9320. [csv R, column 13.]
19. Angelier, J. (1979). Determination of the mean principal
    directions of stresses for a given fault population.
    Tectonophysics 56, T17–T26. [csv phi, column 12.]
20. Simpson, R.W. (1997). Quantifying Anderson's fault types. Journal
    of Geophysical Research 102(B8), 17909–17919. [csv Aphi, col 14.]

Group 4 — material-model provenance (the density column of
Section 2.1):

21. Yeh, T.-Y., Ben-Zion, Y. & Olsen, K.B. (2026). Multi-Scale P and
    S Seismic Velocity Models of California and Western Nevada.
    ESSOAr preprint, doi:10.22541/essoar.15001929/v1.
    [PREPRINT, not yet peer-reviewed — the MUSCAL model behind
    `safs_material_cvm.nc` (CVM Explorer abbr "muscal", UCVM
    component `SCECcode/muscal`, authors Yeh & Ben-Zion). Confirm
    the density-assignment recipe here; replace with the
    peer-reviewed citation when it appears.]
22. Doody, C., Rodgers, A., Afanasiev, M., Boehm, C., Krischer, L.,
    Chiang, A. & Simmons, N. (2023). CANVAS: an adjoint waveform
    tomography model of California and Nevada. Journal of
    Geophysical Research: Solid Earth 128(12), e2023JB027583.
    doi:10.1029/2023JB027583. [Base model MUSCAL builds on;
    verified.]

Group 5 — frictional-strength cap on the differential stress (the
near-surface cap of closure C2, rev 6; Byerlee 1978 = ref 16 and
Zoback & Healy 1992 = ref 2 above also carry this claim). Verified
online (publisher page) 2026-06-15:

23. Sibson, R.H. (1974). Frictional constraints on thrust, wrench and
    normal faults. Nature 249(5457), 542–544. doi:10.1038/249542a0. [Setting-independent mechanics: the
    frictional differential-stress limit per Andersonian regime
    increases ~linearly with depth from ~0 at the free surface —
    supports D_max proportional to Sv_eff.]
24. Brace, W.F. & Kohlstedt, D.L. (1980). Limits on lithospheric stress
    imposed by laboratory experiments. Journal of Geophysical Research
    85(B11), 6248–6252. doi:10.1029/JB085iB11p06248. [Setting-
    independent: the brittle crust's differential stress is bounded by
    Byerlee frictional strength; "a good upper or lower bound to
    observed in-situ stresses ... for pore pressure hydrostatic" — the
    direct justification for capping D at the frictional ceiling.]
25. Townend, J. & Zoback, M.D. (2000). How faulting keeps the crust
    strong. Geology 28(5), 399–402.
    doi:10.1130/0091-7613(2000)028<0399:HFKTCS>2.3.CO;2. [SAF /
    central-southern California: the crust is maintained in frictional-
    failure equilibrium (Byerlee-bounded) — the cap ceiling is where
    the crust actually sits, not merely an upper bound. Companion to
    Townend & Zoback 2004 = ref 5.]

EXCLUDED (Parkfield / creeping central SAF — user directive
2026-06-12; listed so they are not re-introduced by accident):
Hickman & Zoback (2004) GRL 31 L15S12 (SAFOD pilot-hole stress
magnitudes); Lockner et al. (2011) Nature 472 82–85 (SAFOD core gouge
friction).

---

## 8. Decisions

### RESOLVED by user directive 2026-06-15 (rev 5)

1. Closure default: **C2-s** — the spatial Luttrell differential-stress
   field D(x,y) (Section 3 C2-s, Section 4). Flips the rev-3
   recommendation; C1 ratio k = 3.12 retained as `--closure ratio`.
2. Pore pressure: **hydrostatic** regional default (Coyle & Zoback
   1988; Zoback & Healy 1992; Fulton & Saffer 2009). Fault-local
   excess-P_p (Fulton & Saffer profile) deferred to the friction setup.
5. Density: **MUSCAL** lateral-mean rho(z) profile (Section 2.1);
   `constant` kept only for the magnitude check.
6. C2-s data acquisition: **DONE** — the Luttrell-2017 CSM-portal
   export is in `raw_data/luttrell_differential_stress/CSM_data_1781537253009.csv`
   (same 35-column format as the YHSM csv; D = column 16).
7. MUSCAL density recipe: **downgraded to a non-blocking manuscript
   footnote** — use the shipped rho as-is; confirm the Vp->rho relation
   from the Yeh, Ben-Zion & Olsen preprint before publication, not
   before implementation.

### STILL OPEN

3. Axes: Option A (full CSM eigenvectors, keeps San Gorgonio
   transpression) as default vs Option B (strict Andersonian). Plan
   recommends A.
4. Deep taper: off by default; if enabled, z_seis = 15 km (Nazareth &
   Hauksson 2004) or a locally varying map later?
