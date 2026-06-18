# RSSW friction length-scale calculation (pre-implementation report)

Date: 2026-06-10.  Status: IMPLEMENTED — `safs_fault.yaml` (FL=103 blocks:
spatial rs_a/rs_b/rs_srW LuaMap, rs_sl0=0.10, Tnuc_s 24 MPa) and
`parameters.par` (FL=103 + RS constants) written 2026-06-10 per sections 5,
6, 6b; validated against the SeisSol source (RateAndStateInitializer.cpp)
and the tpv104 reference example.  RUNTIME CHECK on first run: the log must
print "RS parameter source (...): f0 1 - muW 1 - b 0" (b from easi); on
older SeisSol builds spatial rs_b silently falls back to the constant
RS_b = 0.0168 in parameters.par, which is wrong away from the hypocenter.

Goal: `safs_seisol_v2_0_0_RSSW` keeps the LSW setup's stress field, CVM
velocity model, mesh, and nucleation machinery, but replaces linear slip
weakening (FL=16) with TPV103/104-style rate-and-state friction with strong
rate weakening (SeisSol FL=103, slip law).  This report determines D_c (the
RS state-evolution distance, SeisSol `rs_sl0` / SCEC `L`) against the L_b
and L_nuc criteria, and flags the downstream parameter decisions the
calculation exposes.

## 1. Inputs

Friction law (SCEC TPV103/104 doc, `tpv104/benchmark_document/SCEC_validation_slip_law.pdf`):

```
tau = f(V, psi) * sigma_n
f(V, psi)  = a * asinh[ V/(2*V0) * exp(psi/a) ]
dpsi/dt    = -(V/L) * [ psi - psi_ss(V) ]            (slip law)
psi_ss(V)  = a * ln[ (2*V0/V) * sinh( f_ss(V)/a ) ]
f_ss(V)    = f_w + ( f_LV(V) - f_w ) / ( 1 + (V/V_w)^8 )^(1/8)
f_LV(V)    = f_0 - (b - a) * ln(V/V0)
psi_ini    = a * ln[ (2*V0/V_ini) * sinh( tau_ini/(a*sigma_n) ) ]
```

Given (user-specified): `f_0 = 0.6`, `V0 = 1e-6 m/s`, `V_init = 1e-12 m/s`
(TPV103 uses V_ini = 1e-16; our 1e-12 LOWERS the direct-effect strength
barrier — see section 6).

Stress (identical to LSW): `safs_initial_stress.yaml`, the
`shmax100_shmin60_sv55_psi69_az23` field (SHmax=100, Shmin=60, Sv=55,
P_p=20 MPa, depth-constant tensor; on-fault values vary per facet).
Fault-wide VW-zone stats (from the projection VTU, 59,781 cells with
depth <= 15 km of 60,658 total):

```
sigma_n_eff [MPa]:  min 37.5   p5 52.1   median 66.3   max 80.0
tau_mag     [MPa]:  median 20.2 (fault-wide), hypocenter 21.2
```

Material (identical to LSW): `safs_material_cvm.nc` (multiscale_statewise
CVM), sampled at fault cell centroids:

```
mu [GPa]:  min 0.4   p5 8.9   median 32.5   max 47.3
```

Hypocenter cell (606971, 3707270, -4965.62; matched at 108 m):

```
depth 5.05 km, sigma_n_eff = 63.9 MPa, mu = 23.4 GPa, nu = 0.277,
tau_0 = 21.2 MPa  (f_apparent = 0.332)
```

Mesh: 500 m fault elements (z0embed _triq, same `safs_mesh.puml.h5`).

## 2. a(z), b(z) from Allison & Dunham (2021), JGR 10.1029/2020JB021394

Their Figure 3a / section 2.4: a and a-b for wet granite (Blanpied et al.
1991, 1995), temperature converted to depth with the 60-km-LAB geotherm.
The paper neglects the shallow VS-to-VW transition (~100 C), so:

```
a - b = -0.004                      for 0 <= depth <= ~15 km   (VW band)
VW -> VS crossing at ~15 km         (350 C for the 60-km LAB geotherm)
a(z): increases ~linearly with depth (linear in T per Rice et al. 2001)
```

The VW->VS crossing at ~15 km conveniently coincides with the LSW deep
barrier (mu_s = 1e6 for z in (-20, -15) km); the RSSW port can realize the
barrier as a TPV103-style velocity-strengthening band instead.

a(z) digitization (already in the repo from the MFEM RS-SRW port,
`friction/rate-and-state/param_a_500m_strongvw.csv`, consistent with
Fig 3a; depth in km, linear interpolation):

```
depth [km]:   0.0       14.892     52.683
a         :   0.006277  0.025508   0.145219
b(z) = a(z) + 0.004  in the VW band
```

At the hypocenter (5.05 km): `a = 0.0128, b = 0.0168, b - a = 0.004`.

Explicit closed-form equations (d = depth in km = -z/1000, z up-positive):

```
a(d):   0 <= d <= 14.892 :  a = 0.0062769 + 0.0012914 * d
        d  >  14.892     :  a = 0.0255077 + 0.0031676 * (d - 14.892)

(a-b)(d):  0 <= d <= 10  :  a - b = -0.004                  (VW band)
           d  >  10      :  a - b = -0.004 + 0.0008*(d - 10)
                            (crosses 0 at 15 km; +0.004 at 20 km)

b(d) = a(d) - (a-b)(d)
```

Checks: a(0) = 0.0063, b(0) = 0.0103; a(5.05) = 0.0128, b(5.05) = 0.0168;
a(15) = 0.0258, b(15) = 0.0258 (neutral); a(20) = 0.0417, b(20) = 0.0377 (VS).

Companion (TPV103 eq 9 analog): ramp V_w from 0.1 to ~1.0 m/s across the
15-20 km band so strong rate weakening is also disabled in the VS barrier.

NOTE: the MFEM RS-SRW SAFS config
(`spatial_friction_rate_state_safs_projected_stress_srw_Dc010_nuc8km_500m.toml`)
uses this same a(z) but a STRENGTHENED `a-b = -0.0134` (its `param_a_minus_b`
CSV) — that was the "invariant lever" used to shrink L_inf for spontaneous
aging-law nucleation.  THIS setup follows Allison exactly (a-b = -0.004),
which fixes the L_nuc/L_b ratio (section 3) and is the crux of the Dc
trade-off.

NOTE: Allison's own d_c = 3.2-10 mm (their Table 1) is NOT transferable:
they size h* = mu*d_c/(sigma_n*(b-a)) = 5 km using the tiny local (b-a)
near their VW->VS crossing at 12 km and resolve it with a much finer grid.
Our floor is the 500 m mesh.

## 3. Formulas and invariants (project standard)

```
L_b   = mu * Dc / (b * sigma_n)                 cohesive-zone scale (mode III;
                                                mode II = L_b/(1-nu), less binding)
L_nuc = mu * Dc / ((b-a) * sigma_n)             Allison's h* (Rice/Ruina style)
L_inf = (2/pi) * mu * b * Dc / ((b-a)^2 * sigma_n)   aging-law critical length
```

Dc-independent invariants with the Allison values (b = 0.0168 at hypocenter):

```
L_nuc / L_b = b/(b-a)              = 4.20
L_inf / L_b = (2/pi)*(b/(b-a))^2   = 11.2
```

So resolving L_b at >= 4 cells (>= 2.0 km) FORCES L_nuc >= 8.4 km and
L_inf >= 22 km.  L_inf cannot fit on this fault — acceptable, because the
slip law has no finite aging-law L_inf instability and nucleation here is
the FORCED overstress patch (TPV103 itself has L_nuc = 26.7 km and
L_inf = 59.5 km against a 3-km patch radius and nucleates by brute-force
overstress).  The binding constraints are therefore:

```
C1 (resolution):  L_b >= 4-5 cells of 500 m   ->  L_b >= 2.0-2.5 km
C2 (nucleation):  L_nuc <= patch DIAMETER (2 x e-fold radius R = 12 km
                  for the LSW patch R = 6000 m), with margin ~1.3 per the
                  proven MFEM aging-law case (8 km patch / L_inf 6.06 km = 1.32)
C3 (containment): L_nuc <= VW depth extent = 15 km
```

## 4. Dc sweep

Hypocenter (mu = 23.4 GPa, sigma_n = 63.9 MPa, b = 0.0168):

```
Dc [m]   L_b [km]  cells@500m   L_nuc [km]   L_inf [km]
0.05      1.09       2.2           4.6          12.2
0.08      1.74       3.5           7.3          19.6
0.10      2.18       4.4           9.1          24.4
0.12      2.61       5.2          11.0          29.3
0.15      3.26       6.5          13.7          36.7
0.20      4.35       8.7          18.3          48.9
```

Fault-wide (all VW cells, per-cell mu/sigma_n/b(z)) — L_b in cells @500 m:

```
Dc [m]   p5     median   frac<3cells   frac<4cells   frac<5cells
0.05     1.28    2.44       82.1%         99.3%        100.0%
0.08     2.05    3.90       11.4%         54.7%         87.4%
0.10     2.56    4.87        6.2%         18.2%         54.7%
0.12     3.08    5.85        4.9%          7.4%         23.1%
0.15     3.84    7.31        3.8%          5.3%          7.4%
0.20     5.13    9.75        2.9%          3.8%          4.9%
```

Where the Dc = 0.10 under-resolved tail lives: the <3-cell set (6.2%) has
median depth 0.5 km (p95 1.6 km), median mu 5.8 GPa — the soft shallow CVM
sediments, where the LSW process zone shrinks identically (Lambda_0 is also
proportional to mu).  Below 3 km depth only 7.3% of VW cells are < 4 cells.

## 5. Recommendation

```
Dc = 0.10 m
```

- C1: hypocenter L_b = 2.18 km = 4.4 cells; fault-wide median 4.9 cells.
  (Slightly under the 5-cell ideal because the CVM hypocenter mu = 23.4 GPa
  is softer than the 32 GPa used in the MFEM sizing; 4.4 cells matches the
  LSW Candidate B level, Lambda_0 = 5.1 cells.)
- C2: L_nuc = 9.1 km vs patch diameter 12 km (R = 6000 m e-fold retained
  from LSW) -> margin 1.31, the same margin class as the proven MFEM
  aging-law case (1.32).
- C3: L_nuc = 9.1 km <= 15 km VW extent.  OK.

Fallback if propagation stalls from under-resolution: Dc = 0.12
(5.2 cells at the hypocenter, <4-cell fraction drops 18.2% -> 7.4%) and
grow the patch e-fold radius R 6000 -> 7000 m to restore the C2 margin
(L_nuc = 11.0 km vs 14 km diameter = 1.27).

## 6. Flags exposed by the calculation (decisions for the config step)

These are NOT part of the Dc determination but fall out of the same numbers;
they need decisions before writing `safs_fault.yaml` / `parameters.par`:

1. **f_w (fully weakened friction) — the project default 0.3 likely kills
   propagation here.**  Ambient dynamic stress drop tau_0 - f_w*sigma_n at
   the hypocenter / median over VW cells deeper than 2 km / fraction of
   those cells with NEGATIVE drop:

   ```
   f_w = 0.20:  8.4 MPa  /  7.1 MPa  /  24% negative
   f_w = 0.25:  5.2 MPa  /  3.8 MPa  /  34% negative   (= LSW mu_d = 0.25 level)
   f_w = 0.30:  2.0 MPa  /  0.5 MPa  /  48% negative
   ```

   Strength-excess ratio using the rate-state direct-effect barrier
   a*sigma_n*ln(V_dyn/V_init) = 22.6 MPa (V_dyn = 1 m/s):

   ```
   S_eff(f_w=0.20) = 2.7    (TPV103 equivalent: 2.8 — known to propagate)
   S_eff(f_w=0.25) = 4.4
   S_eff(f_w=0.30) = 11.4
   ```

   Recommendation: f_w = 0.2 (the TPV103/104 value).  V_w = 0.1 m/s
   (TPV103 and project default).

2. **Nucleation overstress amplitude: the LSW dtau = 10 MPa is NOT enough
   for RS.**  The RS strength the patch must overcome is the direct-effect
   peak ~tau_0 + a*sigma_n*ln(V_dyn/V_init) ~ 21.2 + 22.6 = 43.8 MPa
   (similar to f_0*sigma_n = 38.4 MPa).  TPV103's dimensionless amplitude
   dtau_0/sigma_n = 45/120 = 0.375 maps to ~24 MPa here.  Start with
   dtau = 24 MPa, Gaussian e-fold R = 6000 m, smoothStep ramp t_0 = 1 s
   (same SCEC G(t) machinery already in the LSW yaml/par).

3. **V_init = 1e-12 m/s** (user-specified) enters psi_ini via the formula in
   section 1 (SeisSol computes psi_ini internally from RS_iniSlipRate);
   the choice (vs TPV103's 1e-16) reduces the direct-effect barrier by
   a*sigma_n*ln(1e4) = 7.5 MPa — it helps nucleation and is consistent with
   the MFEM spatial-driver RS initialization.

4. **Deep barrier port**: replace the LSW mu_s = 1e6 band (z in (-20,-15) km)
   with velocity-strengthening a(z) (Allison's natural VW->VS crossing is at
   ~15 km) and/or a TPV103-style delta-a boxcar; below 15 km b = a - |a-b|
   with a-b ramping positive.  Exact ramp shape to be fixed at implementation.

## 6b. Overstress amplitude derivation (accepted params: f_w=0.2, V_w=0.1, Dc=0.10)

Step 1 — the fault starts in frictional equilibrium.  SeisSol sets psi_ini
from eq (11) so that f(V_init, psi_ini)*sigma_n = tau_0 exactly.  Any added
shear stress accelerates V; the question is whether it can push V to seismic
rates (past V_w) before state evolution / stress transfer matter.

Step 2 — frozen-state strength curve.  Accelerating from V_init to V with
psi held fixed (valid while slip << Dc), the asinh law gives, using
asinh(x) ~= ln(2x) for x >> 1:

```
f(V, psi_ini) - f(V_init, psi_ini) = a * ln(V / V_init)
tau_strength(V) = tau_0 + a*sigma_n*ln(V/V_init)        (the direct-effect barrier)
```

Step 3 — barrier amplitude at the hypocenter (a = 0.0128, sigma_n = 63.9 MPa,
V_init = 1e-12 m/s):

```
to reach V_w  = 0.1 m/s:  dtau_min = 0.0128*63.9*ln(1e11) = 20.7 MPa
to reach V    = 1.0 m/s:  dtau_min = 0.0128*63.9*ln(1e12) = 22.6 MPa
```

Step 4 — margin, anchored to TPV103.  TPV103's barrier is
0.010*120*ln(1/1e-16) = 44.2 MPa and it prescribes dtau_0 = 45 MPa — a
1.02x exceedance.  Ours with the same philosophy:

```
dtau_0 = 24 MPa   (= 1.06 x 22.6; also dtau_0/sigma_n = 0.376 vs TPV103's 0.375)
```

Cross-check vs nominal static strength: tau_0 + dtau_0 = 45.2 MPa vs
f_0*sigma_n = 38.4 MPa -> overshoot 0.107*sigma_n; TPV103: (40+45-72)/120
= 0.108*sigma_n.  Identical placement (both problems have tau_0/sigma_n = 1/3).

Step 5 — overstressed footprint with the Gaussian F(r) = exp(-r^2/R^2),
R = 6000 m:

```
r_os(V) = R * sqrt( ln( dtau_0 / (a*sigma_n*ln(V/V_init)) ) )
r_os(V_w = 0.1 m/s) = 6000*sqrt(ln(24/20.7)) = 2.31 km
r_os(1.0 m/s)       = 6000*sqrt(ln(24/22.6)) = 1.47 km
```

TPV103's equivalent driven-to-seismic core is only ~0.4 km (its compact
F(r) falls below 44.2/45 = 0.982 at r = 0.133*R = 400 m) and propagates
fine — once the core goes seismic and weakens toward f_w, the crack-tip
stress concentration takes over; sustained propagation is governed by the
ambient S_eff = 2.7 (section 6.1), not by the footprint.  Our 1.5 km core
is ~4x TPV103's.

Conservatisms not credited: state evolution during the 1 s ramp (weakens,
lowers the barrier) and quasi-static stress transfer from patch slip
(concentrates stress at the patch edge).  Both help; 24 MPa is mildly
conservative, like TPV103's 45.

### 6b.1 The LSW onset check rewritten for rate-and-state (yaml-comment-ready)

LSW original (safs_fault.yaml history): patch stress should exceed static
strength where the Gaussian is over the strength excess SE = mu_s*sn - tau_0:
r_os = R*sqrt(ln(dtau/SE)) = 6000*sqrt(ln(10/6.24)) = 4122 m > a_c(Uenishi-
Rice, LSW) = 2289 m -> onset ratio 1.80.

RS-SRW version (hypocenter mu=23.4 GPa, nu=0.277, sn_eff=63.9 MPa,
tau_0=21.2 MPa; f0=0.6, V0=1e-6, a=0.0128, b=0.0168, Dc=0.10, f_w=0.2,
V_w=0.1, V_init=1e-12; dtau_nuc=24 MPa, R=6000 m e-fold):

```
RS strength excess (direct-effect barrier to V_dyn = 1 m/s; replaces
LSW's SE = mu_s*sn - tau_0 — RS has no static mu_s, strength is rate-
dependent and the patch must out-climb the direct effect at frozen state):
  SE = a*sn*ln(V_dyn/V_init) = 0.0128*63.9*ln(1e12) = 22.6 MPa

overstressed footprint (same Gaussian algebra as LSW):
  dtau_nuc*exp(-r_os^2/R^2) = SE
  -> r_os = R*sqrt(ln(dtau_nuc/SE)) = 6000*sqrt(ln(24/22.6)) = 1471 m
  (to V_w = 0.1 m/s instead: SE = 20.7 MPa -> r_os = 2301 m)

critical half-length — NOT the quasi-static RS length (L_nuc/2 = 4574 m;
r_os can never beat that, and TPV103 itself sits at r_os = 396 m vs
L_nuc/2 = 13.3 km).  Inside r_os the fault is driven past V_w, so the
crack that must go supercritical is the SRW-WEAKENED one; Uenishi-Rice
with the slip-law breakdown slope W = (f_peak - f_ss(V_dyn))*sn/Dc:
  f_peak       = (tau_0 + SE)/sn = 0.685
  f_ss(1 m/s)  = f_w + (f_LV(1) - f_w)/(1+(1/V_w)^8)^(1/8)
               = 0.2 + (0.5447-0.2)/10.0 = 0.2345
  a_c = 0.579*mu*Dc/((1-nu)*(f_peak - f_ss)*sn)
      = 0.579*23.4e9*0.10/(0.723*0.451*63.9e6) = 65 m

onset ratio r_os/a_c = 1471/65 = 22.6
calibration: the same check applied to TPV103's published design gives
r_os = 396 m, a_c = 176 m, ratio 2.2 (LSW criterion here was 1.80) —
so >~ 2 is the empirically-blessed regime; 22.6 is ample.
```

Resolution caveat exposed by the same numbers: the SRW dynamic breakdown
zone Lambda_0_srw = (9pi/32)*mu*(2*Dc)/((f_peak-f_ss)*sn) ~= 143 m, i.e.
SUB-CELL at h = 500 m (TPV103's is 404 m, resolved by its 100-200 m
meshes).  The L_b >= 4-5 cell criterion (section 5) guards the quasi-
static state-evolution scale; the dynamic SRW front will be numerically
spread over ~2-3 cells.  Same regime the MFEM RS-SRW 500 m sibling
accepted (its config marks the run "propagation/stability smoke test");
sharpening it means h ~ 100-150 m near-fault or a larger Dc (which
trades directly against L_nuc and the patch margin).

## 7. Reproduction

All numbers from a Python session over:
- `stress/results/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax100_shmin60_sv55_psi69_az23/*_fault_stress.vtu` (cell data, 60,658 tris)
- `seisol_quakeworx/safs_seisol_v2_0_0_RSSW/safs_material_cvm.nc` (trilinear mu at centroids)
- a(z): `friction/rate-and-state/param_a_500m_strongvw.csv`, b = a + 0.004, VW = depth <= 15 km

Hypocenter matched to nearest fault cell centroid (108 m).

## 8. TEST2 run result + arrest analysis (2026-06-11) -> f_w = 0.1

TEST2 (Expanse, 2 nodes, EndTime 100 s, output every 5 s, SRs/SRd/Sls/Sld
only) NUCLEATED AND RAN: M0 rate peaked 1.2e19 N m/s at t~5 s, event over
by t~10 s at Mw 7.21.  Rupture extent (slip > 0.1 m): hit the SW fault end
(-18.4 km) and ARRESTED at +43 km NE of a 275-km strand.  The DR surface in
the puml mesh spans depth 0 to -16.6 km (NOT the full -41.6 km of the
stress-projection fault surface).

Along-strike stress profile (VW core 2-15 km depth, from the projection
VTU) explains the arrest:

```
along [km]    mu_app(med)   drop tau0-0.2*sn [MPa]
-20..+40      0.33 -> 0.29   8.5 -> 5.9     <- ruptured
+40..+70      0.245-0.254    3.3 - 3.8      <- BARRIER (arrest 3 km in)
+80..+170     0.34 -> 0.40   8.8 - 11.4     <- strongest stretch
+180..+190    0.16 - 0.23    ~0 to -3       <- field dying
+200..+275    0.02 - 0.04    -13 to -14     <- tractionless, unbreakable
```

G_c ~ (f_peak - f_ss)*sn*Dc ~ 2.9 MJ/m2 in the corridor vs ambient drop
3.4 MPa -> energy deficit -> arrest.  Nucleation NOT implicated (worked;
patch/amplitude unchanged).

DECISION (user, 2026-06-11): f_w = 0.1 (RS_muW 0.2 -> 0.1).
- corridor drop 3.4 -> 10.5 MPa; hypocenter drop 14.8 MPa, S_eff = 1.53
  (< 1.77 Burridge-Andrews: SUPERSHEAR likely in the strong 80-170 km
  stretch — expected feature, not a bug).
- positive-drop coverage now extends to ~190 km; beyond that the field is
  tractionless and no friction setting can carry rupture (the only lever
  there is the regional stress itself, az 23 / SHmax 100 / Shmin 60).
- full -18..+180 km rupture ~ 200 km ~ Mw 7.8.

Known side effects (accepted): slip grows further (TEST2 already peaked
~50 m, concentrated in soft CVM sediments mu 0.4-6 GPa, dip-heavy rake
~-77 deg from the stress projection; broken zone relaxes at ~1-2 m/s for
tens of s, partly fed by the PERMANENT +24 MPa patch).  Slip-magnitude
realism is a separate calibration if QuakeWorx needs it.
