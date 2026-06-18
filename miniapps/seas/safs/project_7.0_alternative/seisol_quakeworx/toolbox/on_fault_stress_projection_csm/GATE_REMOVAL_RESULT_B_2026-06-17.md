# Gate-removal result — recipe B (end-to-end, baked nc)

Sampled the baked nc at 60658 PUML fault facet centroids (as SeisSol does) and computed S = (mu_s - mu_app)/(mu_app - mu_d).

- BEFORE: `safs_stress_csm.nc` with mu_s=0.47 (production V3).
- AFTER : `safs_stress_csm_B.nc` with mu_s=0.318 (recipe B: lower k + fault-local overpressure).
- mu_d=0.1, S_TARGET=1.7, seismogenic z in (-12.0, -3.0) km.

## Per-gate seismogenic S < 1.7 (before -> after)
| gate | n_seis | mu_app med B->A | pure-barriers B->A | %S<1.7 B->A | residual S>=1.7 facets (after) |
|---|---|---|---|---|---|
| SanGorgonio | 4360 | 0.191->0.242 | 158->4 | 18%->98% | 69 |
| g181 | 615 | 0.176->0.239 | 0->0 | 0%->97% | 16 |
| g221 | 2766 | 0.198->0.242 | 0->0 | 2%->100% | 13 |
| g266 | 0 | - | - | - | (not on PUML mesh) |

## Global checks on the AFTER field
- Subcriticality: max mu_app (active z in (-15.0, -0.3) km) = 0.2971  <  mu_s=0.318  -> PASS (no t=0 pre-slip).
- Tension: min sigma_n_eff (all facets) = 0.20 MPa; min in seismogenic band = 18.38 MPa  -> PASS (>0).
- Overall: ALL TARGETED GATES REMOVED (>=95% S<1.7).

NOTE: S<1.7 is the analytic propagation proxy (Das & Aki). A SeisSol dynamic-rupture run on Frontera is the ultimate confirmation that the rupture crosses; this static check is necessary, not sufficient.

## Nucleation check (hypocenter, recipe B)
The hypocenter (s=0, SE) is outside the gate bands, so the overpressure does not
touch it; only k=1.8 changes its stress. The existing 20 MPa Gaussian patch
(R=6 km) still nucleates:
- sigma_n_eff 71.5 MPa, tau_0 15.4 MPa, mu_app 0.215; tau_s=mu_s*sigma_n=22.7 MPa,
  strength excess SE=7.4 MPa.
- r_os = R*sqrt(ln(20/7.4)) = 5998 m;  Uenishi-Rice Lnuc = 6017 m;  Lnuc/2 = 3009 m.
- r_os 5998 > Lnuc/2 3009  ->  onset ratio 1.99 > 1  ->  NUCLEATES (production
  was 2.33; tighter because mu_s-mu_d shrank 0.37->0.218, but still > 1).

## Residual facets (not fully removed)
A small number of deeply-clamped barrier cores cannot reach the floor within the
lambda<=0.9 sub-lithostatic cap (SanGorgonio 69 / g181 16 / g221 13, all < 2.5%
of the band; SanGorgonio retains 4 pure barriers mu_app<mu_d). They are isolated
relative to the process zone Lb (~few km here), so a dynamic rupture is expected
to cross them; the SeisSol run will confirm.

## Production staging (NOT yet promoted)
The production files are UNTOUCHED. Recipe B is staged alongside them for
inspection/diff:
- `safs_seisol_v3_0_0_LSW/safs_stress_csm_B.nc`        (baked gate-removed nc, 79.6 MB)
- `safs_seisol_v3_0_0_LSW/safs_initial_stress_B.yaml`  (ASAGI -> the _B nc)
- `safs_seisol_v3_0_0_LSW/safs_fault_B.yaml`           (mu_s=0.318, includes _B stress)
- `toolbox/.../overpressure_B_2026-06-17.npz`          (the (s,z) overpressure field)

To run V3 with recipe B, point `parameters.par` at `safs_fault_B.yaml` (or, to
promote, rename the _B files over the originals after backing them up). To try
recipe C instead, rebuild the field with `make_overpressure_field.py --k-ratio
1.5 --recipe C --target-mu 0.18` and re-bake with `--k-ratio 1.5` + mu_s 0.225.
