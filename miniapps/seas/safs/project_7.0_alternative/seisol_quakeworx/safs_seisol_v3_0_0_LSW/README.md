# safs_seisol_v3_0_0_LSW — CSM depth-dependent stress + constant LSW friction

SeisSol single-event dynamic rupture (FL=16 linear slip weakening) on the SAFS
multi-strand fault. V3 = V2's case structure, but the initial fault stress is
the **depth-dependent CSM C1 k=2.39 projection** (instead of V2's single
constant regional tensor), and the friction is **constant** mu_s / mu_d / d_c.

Full design + analysis: `PLAN.md` (+ `PLAN.pdf`). Diagnostics:
`csm_S_ratio_const_mu_strike_depth.png`, `csm_S_ratio_const_mu_fault.vtu`.

## What's in this folder
| file | role |
|---|---|
| `parameters.par` | SeisSol namelists (FL=16, godunov, LTS, fault+energy output, EndTime 100 s) |
| `safs_initial_stress.yaml` | initial stress — `!ASAGI safs_stress_csm.nc` (+ ConstantMap fallback) |
| `safs_stress_csm.nc` | the spatially varying stress (compound `data{s_xx..s_xz}`, comp-NEGATIVE Pa, 1000 m x/y, 250 m z, 80 MB) |
| `safs_fault.yaml` | friction: **constant** mu_s=0.47 (LuaMap, +deep barrier) / mu_d=0.10 / d_c=2.5 / cohesion=0; Tnuc_s Gaussian nucleation |
| `safs_material_cvm.yaml`, `safs_material_cvm.nc` | MUSCAL CVM material (symlink → V2; unchanged) |
| `safs_mesh.puml.h5` | PUML mesh (symlink → V2; unchanged) |

Regenerate the stress nc:
```
conda run -n pythonenv python3 \
  ../toolbox/on_fault_stress_projection_csm/csm_stress_to_asagi.py \
  --csm-csv ../raw_data/yang_and_hauksson_orientation/CSM_data_1781290479811.csv \
  --material-nc ../safs_seisol_v2_0_0_RSSRW/safs_material_cvm.nc \
  --mesh ../safs_seisol_v2_0_0_RSSRW/safs_mesh.puml.h5 \
  --verify-vtu ../toolbox/on_fault_stress_projection_csm/csm_yhsm2013_stress_on_safs_mesh_fault_stress.vtu \
  --out safs_stress_csm.nc --dx 1000 --dz 250
```

## Friction values and why
- **mu_s = 0.47** — must exceed the fault-wide max mu_app = 0.449 to avoid
  spontaneous t=0 slip (the nc-interpolated max is 0.448, so 0.47 leaves a
  0.021 margin). Set in the `mu_s` LuaMap, with the deep barrier mu_s=1e6 for
  z in (−20000, −15000) m.
- **mu_d = 0.10** — SAF weak-fault dynamic friction (heat-flow / low apparent
  friction). ~3% of facets have mu_app < 0.10 (gate cores + shallow patches);
  there mu_d ≥ mu_app so they never weaken (pure barriers).
- **d_c = 2.5 m** — process zone Lb ≈ 2.3 km (~4.5 elements) at the hypocenter;
  marginal (~1 element) only near the deep barrier.
- nucleation: Tnuc_s = 20 MPa Gaussian, R = 6000 m, at (606971, 3707270,
  −4965.62), ramped over [0,1] s.

## EXPECTED OUTCOME — read before interpreting the run
With the CSM stress, mu_app varies 0.014–0.449 along strike, and a **constant**
mu_s pinned above 0.449 makes the low-mu_app **San Gorgonio gate** (s ≈ 40–60
km, mu_app ≈ 0.14, SHmax ≈ 78° to the fault strike — confirmed by both YHSM and
Luttrell) super-critical: S = (mu_s−mu_app)/(mu_app−mu_d) ≈ 14 there, far above
the ~1.77 propagation limit, and it is a full-depth wall ~24 km long even at the
free surface. **The SE-nucleated rupture is therefore expected to ARREST at the
San Gorgonio gate, not reach the NW end.** This is the faithful, physically
realistic result of constant friction on the real CSM stress (the SAF segments
at San Gorgonio); it is *not* a full SE→NW rupture. Full propagation under
constant friction would require a fault-favorable (near-uniform mu_app) stress
that departs from the CSM orientation — see PLAN.md §4.1 and the trilemma.

## Run notes (cluster)
- Requires SeisSol built with **ASAGI** (the QuakeWorx app has it).
- The three symlinked files (`safs_mesh.puml.h5`, `safs_material_cvm.nc`,
  `safs_material_cvm.yaml`) must be **materialized as real files** before
  transferring to the cluster (`cp -L`), or transfer V2's copies alongside.
- `mkdir output` (or let the job create it) for the `output/safs` prefix.
- Mesh/material/numerics and the run command are identical to V2 — reuse the
  V2 sbatch with this folder as the case dir.

## Validation gates
- [x] Stress nc self-check (LOCAL, 2026-06-15): trilinear samples reproduce the
      C1 projection VTU sigma_n_eff / tau to median 1.4e-3 / 2.0e-3 (p95 2–3%);
      nc-interpolated max mu_app 0.448 < mu_s 0.47 → no t=0 pre-slip.
- [ ] (cluster) t=0 fault output Ts0/Td0/Pn0 matches the projection VTU
      per-facet (OutputMask group 5 is on).
- [ ] (cluster) nucleation initiates at the hypocenter; cohesive zone resolved.
- [ ] (cluster) rupture front (RT, group 10) propagates SE and arrests at the
      San Gorgonio gate (s ≈ 40–60 km) — the expected stopping point.
- [ ] (cluster) reference-vector sign check (XRef/YRef/ZRef) on the curved fault
      (V2 REVIEW.md R-001): magnitudes are convention-independent; if T_s/T_d
      signs look flipped, fix via the reference vector, not the stress.
