# MFEM SAFS → SeisSol porting map

Source problem: `seas-mfem-safs/miniapps/seas` SAFS LSW dynamic-rupture run, job
`spatial_dyn_slipweakening_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch`,
config `…/safs/project_7.0_alternative/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml`.

Target: a SeisSol single-event dynamic-rupture case (FL=16 linear slip weakening),
expressed as `parameters.par` + `safs_fault.yaml` + `safs_initial_stress.yaml` +
`safs_material.yaml`, in the updated **tpv13 / QuakeWorx** file format.

---

## 0. The physics, in one paragraph

A real San-Andreas-Fault segment (SAF-ALT6), meshed in **UTM 11N** coordinates
(x=East, y=North, z=Up; z<0 = depth). Constant elastic medium (λ=μ=32 GPa,
ρ=2670 → Vp≈5.4 km/s, Vs≈3.46 km/s, ν=0.25). **Linear slip-weakening** friction
(μ_s=0.85, μ_d=0.30, D_c=2.0 m, cohesion=0). Background stress = a regional
**Hickman & Zoback constant Cauchy tensor** (compression POSITIVE), effective
normal stress via **pore pressure P_p=16 MPa**. Rupture is nucleated by a
**time-ramped Gaussian over-stress** (Δτ_strike=20 MPa, radius 4 km, ramp T_nuc=1 s)
at hypocenter (606971, 3707270, −4965.62). A **deep barrier** locks z∈[−20,−15] km.

---

## 1. Convention differences (the high-risk part)

| Aspect | MFEM SAFS | SeisSol | Action |
|---|---|---|---|
| Stress sign | compression **positive** | compression **negative** | **negate the whole tensor** |
| Axes | E,N,U (UTM 11N) | x=E, y=N, z=U | identical → component-wise 1:1 |
| Strike traction sign | right-lateral positive | right-lateral positive (`docs/dynamic-rupture.rst:31`) | same (validate) |
| Pore pressure | `σ_n_eff = σ_n_total − P_p` applied by driver | **not auto-applied** by LSW | bake P_p into the stress diagonal |
| σ_n strength floor (10 MPa) | MFEM stabilizer (`max(σ_n,floor)`) | **no equivalent** | drop (see §6); not a physical knob |
| ADER order | `ader_order = 2` (runtime) | **compile-time** `CONVERGENCE_ORDER` | build flag, not in .par |
| Flux | `mixed_flux=none` (pure upwind) | Godunov is the default upwind flux | `numflux='godunov'` |

SeisSol LSW strength (`src/DynamicRupture/FrictionLaws/CpuImpl/LinearSlipWeakening.h:153-160`):
`strength = −cohesion − μ·min(σ_n_total, 0)` with σ_n negative in compression.

---

## 2. Stress tensor port + projection (the numbers)

Effective Cauchy tensor that MFEM projects (compression POSITIVE, E-N-U), with
P_p already removed from the diagonal:

```
σ_total (MPa, comp+):  xx=50.566  yy=111.434  zz=45.0  xy=−9.889  yz=0  xz=0
P_p = 16 MPa
σ_eff   (MPa, comp+):  xx=34.566  yy=95.434   zz=29.0  xy=−9.889  yz=0  xz=0
```

SeisSol values = **negate** (comp− convention), placed in `safs_initial_stress.yaml`:

| SeisSol field | Value (Pa) | Note |
|---|---|---|
| `s_xx` | −3.4566191478555e7 | 34.566 MPa eff., negated |
| `s_yy` | −9.5433808521445e7 | 95.434 MPa eff., negated |
| `s_zz` | −2.9000e7 | 29.0 MPa eff., negated |
| `s_xy` | +9.888543819998e6 | shear, sign-flipped from −9.889 (comp+) |
| `s_yz` | 0 | |
| `s_xz` | 0 | |

**SeisSol projects this Cartesian tensor onto the fault — confirmed in source.**
`BaseDRInitializer::rotateStressToFaultCS` (`BaseDRInitializer.cpp:127,249`) applies
`inverseSymmetricTensor2RotationMatrix(fault.normal, tangent1, tangent2)` per facet, i.e.
`σ_n = n·(σ·n)`, `τ_strike/dip = t·(σ·n)`. This is the **identical** projection MFEM does
in `FaultGeometry::ComputeSAFSParams` (`fault_geometry_safs_templated.inl:102-108`:
`sigma_n_total = n·Sn; tau1 = t1·Sn; tau2 = t2·Sn; sigma_n_eff = sigma_n_total − P_p`, with
`Sn = σ·n`). So a constant tensor is the correct input; **do not pre-project**, and on the
curved fault the per-facet tractions vary automatically.

Why baking P_p into the diagonal is exact: P_p is isotropic, so projecting
`(σ_total − P_p·I)` onto the fault normal gives `σ_n_total − P_p` (= σ_n_eff), and
the tangential projection of `P_p·I` is zero → shear unchanged.

SeisSol alternatives (also supported, not needed for a constant tensor): fault-local
tractions `T_n/T_s/T_d`, or a gridded field via easi `!ASAGI` (HDF5).

**Validation gate** (after first run, fault output Ts0/Td0/Pn0 at the hypocenter DOF,
comp+ equivalents): MFEM reports σ_n_eff≈49.27 MPa, τ_strike0≈26.32, τ_dip0≈−12.65,
|τ_pre|≈29.2. Magnitudes are convention-independent and should match. If |Pn0| matches but
T_s/T_d signs are flipped, fix via the reference vector (§4), **not** by editing stresses.

---

## 3. Friction → fault YAML (FL=16)

Required FL=16 fields (`LinearSlipWeakeningInitializer.cpp:58-66`): `mu_s`, `mu_d`,
`d_c`, `cohesion`; optional `forced_rupture_time`.

| MFEM | value | SeisSol field | value |
|---|---|---|---|
| `mu_s_default` | 0.85 | `mu_s` | 0.85 |
| `mu_d_default` | 0.30 | `mu_d` | 0.30 |
| `d_c_default` | 2.0 | `d_c` | 2.0 |
| `cohesion_default` | 0.0 | `cohesion` | 0.0 (SeisSol cohesion is ≤0; 0 is fine) |
| deep barrier z∈[−20k,−15k], mu_s=1e6 | lock | `mu_s` via LuaMap = 1.0e6 in band | locks (strength → huge) |

`mu_s` is a `!LuaMap` returning 1.0e6 for −20000<z<−15000 else 0.85.

---

## 4. Reference vector & fault frame (XRef/YRef/ZRef, refPointMethod)

`docs/dynamic-rupture.rst:38-60`: the reference point/direction orients the fault normal
and hence strike `s=(n_y,−n_x,0)` / dip. We use `refPointMethod=1`, `(0,−1,0)` (TPV33-style).
**This is the one item to validate at runtime** — for the curved SAF the projected strike/dip
SIGNS depend on this basis (the magnitudes do not). See REVIEW.md R-001.

---

## 5. Nucleation → Tnuc_* + s_0/t_0

MFEM `gradual_overstress`: Gaussian Δτ_strike=20 MPa, radii 4 km, smoothStep ramp over
T_nuc=1 s, center (606971, 3707270, −4965.62).

SeisSol analog (`docs/dynamic-rupture.rst:368-393`; applied as `initialStress +=
Tnuc·smoothStepIncrement(t−s0, dt, t0)` over `[s0, s0+t0]`, `FrictionSolverCommon.h:433-435`):

| MFEM | SeisSol |
|---|---|
| Δτ_strike = 20 MPa, Gaussian | `Tnuc_s` = 20e6·exp(−r²/4000²) via LuaMap |
| Δτ_dip = 0 | `Tnuc_d` = 0 |
| (normal unperturbed) | `Tnuc_n` = 0 |
| ramp start = 0 | `s_0 = 0.0` (`&DynamicRupture`) |
| T_nuc = 1.0 s | `t_0 = 1.0` (`&DynamicRupture`) |

`forced_rupture_time = 1e10` is **included** (tpv13 style) and inert: with 1e10 the LSW
forced-rupture term `f2 = clamp((t−1e10)/t_0,0,1) = 0` forever (`LinearSlipWeakening.h:197-208`),
and it does not touch the independent `Tnuc_s` smoothStep path. Gaussian uses 3-D distance
`r²=(x−606971)²+(y−3707270)²+(z+4965.62)²` (≈ in-plane for on-fault DOFs, equal radii).

---

## 6. The "normalcap" (σ_n strength floor) — NOT ported

MFEM's `sigma_n_strength_floor_pa = 10 MPa` (the "normalcap") is a **numerical stabilizer**
for an MFEM-specific tensile free-slip blow-up, not a physical parameter. SeisSol's ADER-DG LSW
(`min(σ_n,0)`, SCEC-verified) has no such failure mode and **no** equivalent knob — do not
replicate. SeisSol's analogous lever is `etaDamp` (`&DynamicRupture`, <1.0). MFEM's
`min_sigma_n_pa = 1 MPa` clamp (`fault_geometry_safs_templated.inl:110`) is likewise
MFEM-internal; σ_n on the fault is well above 1 MPa except at shallow tips.

---

## 7. Material → material YAML

Constant medium → `!ConstantMap`: `rho=2670`, `mu=3.2e10`, `lambda=3.2e10`. No off-fault
plasticity (`Plasticity=0`), so unlike tpv13 the material file omits `plastCo`/`bulkFriction`
and the volume-stress `!Include`.

**CVM variant (this folder, `safs_seisol_cvm`)**: the constant `!ConstantMap` above is the
verified baseline (`safs_material.yaml`, kept for the A/B switch); the active material is
the **multiscale_statewise CVM** (= MFEM `use_sidecar=true`) via easi `!ASAGI` in
`safs_material_cvm.yaml`, fed by `safs_material_cvm.nc` generated with
`convert_cvm_to_asagi.py` from the MFEM sidecar. See README_port.md, "Heterogeneous
material (CVM)" for the workflow, validation gates, and caveats. (Side note: the
constant medium's Vp is √((λ+2μ)/ρ) ≈ **5996 m/s** — earlier comments saying ≈5.39 km/s
were arithmetically wrong; Vs ≈ 3462 m/s is correct.)

---

## 8. Numerics / discretization / output → parameters.par

| MFEM | SeisSol `.par` | Note |
|---|---|---|
| `cfl=0.5` | `&Discretization CFL=0.5` | direct |
| `ader_order=2` | build-time `CONVERGENCE_ORDER` | not a .par key; ≥3 for fidelity |
| `mixed_flux=none` (upwind) | `&equations numflux='godunov'`, `numfluxnearfault='godunov'` | upwind |
| `use_pml=false` | absorbing BC via mesh tag 5 | no PML in SeisSol |
| `tfinal=100s` | `&AbortCriteria EndTime=100.0` | ABC reflections leak after ~7s |
| `dt_max=0.1s` | `&Discretization FixTimeStep=0.1` | CFL-limited anyway |
| LTS | `&Discretization ClusteredLTS=2` | like tpv13 |
| paraview fault dt 1.0 s | `&Elementwise printtimeinterval_sec=1.0` | 11-field `OutputMask` |
| volume/energy output | `&Output Format=6`/`iOutputMask`/`EnergyOutput=1`/`Checkpoint=0` | tpv13 style |

---

## 9. Mesh — the real blocker (separate from param files)

MFEM mesh: gmsh `.msh` (UTM coords), physical tags **fault=101, top=102, bottom=103,
sides=104, rock=1**. SeisSol needs **PUML** `.puml.h5` (via **pumgen**), boundary codes
**1=free surface, 3=dynamic rupture, 5=absorbing** (`fault-tagging.rst:12`):

| MFEM tag | SeisSol code |
|---|---|
| 101 fault | **3** |
| 102 top | **1** |
| 103 bottom | **5** |
| 104 sides | **5** |
| 1 rock | volume group |

Retag `.geo`/`.msh` to {1,3,5}, then `pumgen <prefix>.msh -s msh2`, set
`MeshFile='safs_mesh.puml.h5'`. **pumgen is not installed** (prior session) — out of scope
for the param-file deliverable. See [[tpv33-mesh-workflow-deps]].

---

## 10. Generated files (in `safs_seissol/`)

- `parameters.par` — namelists per §8 (tpv13 format).
- `safs_fault.yaml` — friction + barrier (Lua) + Tnuc_s (Lua) + forced_rupture_time=1e10;
  `[s_xx..s_xz]: !Include safs_initial_stress.yaml`.
- `safs_initial_stress.yaml` — constant effective stress tensor (§2), SeisSol projects it.
- `safs_material.yaml` — ConstantMap ρ/μ/λ (§7).
- `README_port.md` — run guide + validation gates + mesh prerequisite.

## Open questions / to validate at runtime
1. Sign of T_s/T_d at the hypocenter vs MFEM (magnitudes match; signs depend on the reference vector). R-001.
2. SeisSol build `CONVERGENCE_ORDER` (2 to match MFEM, or ≥3).
3. Mesh conversion + boundary retag (pumgen) — prerequisite.
