# Code Review: SAFS → SeisSol parameter-file port (tpv13-aligned)

## Review Scope
- Plan: `MFEM_to_SeisSol_mapping.md`
- Files reviewed: `parameters.par`, `safs_fault.yaml`, `safs_initial_stress.yaml`,
  `safs_material.yaml`, `README_port.md`
- Domain context: SeisSol source `/Users/chunhuizhao/projects/SeisSol`
  (`src/Initializer/Parameters/*`, `src/DynamicRupture/Initializer/*`,
  `src/DynamicRupture/FrictionLaws/*`, `docs/*.rst`); the updated tpv13 reference
  in `seisol_quakeworx/tpv13/` and `docs/tpv13.rst`; MFEM SAFS source.

## Verified against source (NOT bugs)
- **Friction names** `mu_s, mu_d, d_c, cohesion` exact (`LinearSlipWeakeningInitializer.cpp:58-61`);
  Cartesian `s_xx..s_xz` exact (`BaseDRInitializer.cpp:343`); nucleation `Tnuc_n/Tnuc_s/Tnuc_d`
  exact (`:330`).
- **Stress projection**: SeisSol projects the global Cartesian tensor onto the fault per facet
  (`BaseDRInitializer::rotateStressToFaultCS`, `inverseSymmetricTensor2RotationMatrix`,
  `:127,249`) — the same `n·σ·n` / `t·σ·n` projection as MFEM `ComputeSAFSParams`
  (`fault_geometry_safs_templated.inl:102-108`). The constant-tensor input is correct.
- **Nucleation** ramp = `smoothStepIncrement(t−s0, dt, t0)` over `[s0,s0+t0]`
  (`FrictionSolverCommon.h:433-435`) — faithful to MFEM smoothStep; `nucleationCount` default 1
  (`DRParameters.cpp:78`) so the single `Tnuc_*` patch IS read.
- **`forced_rupture_time = 1e10`** (now included, tpv13 style) is inert: `f2 =
  clamp((t−1e10)/t_0,0,1)=0` forever (`LinearSlipWeakening.h:197-208`), independent of the
  `Tnuc_s` path. Safe.
- **`numflux`/`numfluxnearfault`** read from `&equations` (`ModelParameters.cpp:24,106,114`).
- **Format alignment** matches tpv13: separated `!Include` stress, `.puml.h5` MeshFile,
  11-field `OutputMask`, empty `&Pickpoint`, `Format=6`/`iPlasticityMask`/`Checkpoint=0`.
- **Stress arithmetic** (effective = total − P_p, negated) correct; physics self-consistent
  (μ_app≈0.59 < μ_s ⇒ sub-yield background; +20 MPa nucleation ⇒ ruptures).
- All YAML parses; unknown `.par` keys only warn (`ParameterReader.cpp:76`).

## Findings

### [R-001] [MODERATE] [POSSIBLE] [parameters.par:&DynamicRupture] — Reference vector on a curved fault

**Category:** ASSUMPTION

**Description:** `XRef/YRef/ZRef = (0,−1,0)`, `refPointMethod = 1` is the planar-TPV33 choice.
On the **curved** SAF the strike basis `s=(n_y,−n_x,0)` can flip sign where the normal's
y-component changes sign → reversed `T_s`/`T_d` (and slip direction) on part of the fault.
The projection MAGNITUDES (σ_n, |τ|) are unaffected; only the SIGNS are at risk.

**Suggested fix:** prefer a reference *point* on one side of the whole fault, then validate:
```diff
- XRef = 0.0
- YRef = -1.0
- ZRef = 0.0
- refPointMethod = 1
+ XRef = 1.0e7      ! point far EAST of the whole fault (UTM 11N); tune to the mesh
+ YRef = 3.707e6
+ ZRef = 0.0
+ refPointMethod = 0
```

**Test case (runtime, needs the mesh):**
```
# ParaView on output/safs-fault.xdmf:
#   Pn0 < 0 everywhere (compression); Ts0 single-signed along strike (no checkerboard flip).
# A sign-flip stripe along strike => frame inconsistency => change the reference.
```

---

### [R-002] [LOW] [parameters.par:&Output] — Format=6 volume output is heavy on the SAFS mesh

**Category:** ASSUMPTION

**Description:** `Format = 6` (hdf5 wavefield) matches tpv13, but on the ~3.7 M-tet SAFS mesh
~100 volume snapshots is large. Intentional (style parity), flagged for the operator.

**Suggested fix:** set `Format = 10` for a first smoke (volume output off); fault + surface +
energy output remain.

**Test case:** N/A (output-cost preference).

---

### [R-003] [LOW] [parameters.par:&AbortCriteria] — EndTime 100 s outlives absorbing-BC validity

**Category:** ASSUMPTION

**Description:** Box absorbing BCs (code 5) leak reflections after ~7 s (min_box_dim/cp), as the
MFEM header notes. Documented in README/.par.

**Suggested fix:** `EndTime = 15.0` for a clean first window; raise for the arrest test.

**Test case:** N/A.

---

### [R-004] [LOW] [safs_fault.yaml:Tnuc_s] — Gaussian uses 3-D distance vs MFEM in-plane

**Category:** ASSUMPTION

**Description:** Equal radii (4000 m) + on-fault DOFs ⇒ 3-D distance ≈ in-plane distance;
difference negligible. Documented.

**Suggested fix:** none.

**Test case:** N/A.

---

### [R-005] [LOW] [README/&Output] — `output/` directory must pre-exist

**Category:** EDGE_CASE

**Description:** `OutputFile = 'output/safs'`; SeisSol aborts at write time if `output/` is
missing. README says `mkdir -p output`.

**Suggested fix:** keep the `mkdir -p output` step prominent.

**Test case:** N/A.

---

## Summary
- Critical issues: 0 (the global force-rupture candidate from shared `t_0` is verified absent:
  `f2=0` with `forced_rupture_time=1e10`).
- Moderate issues: 1 (R-001 reference vector — runtime validation; primary action item).
- Low issues: 4.
- Plan compliance: **FULL**; format matches the updated tpv13 reference.
- Verdict: **PASS WITH FIXES** — files are physically and syntactically correct against the
  SeisSol source. R-001 needs the converted mesh to validate. R-002–R-005 are caveats/preferences.

## Unreviewed Areas
- **Mesh** (`safs_mesh.puml.h5`): not produced (pumgen not installed). Boundary retag
  101→3 / 102→1 / 103,104→5 and right-handedness unverifiable until conversion; R-001 also
  depends on the final geometry.
- **SeisSol build** (`CONVERGENCE_ORDER`, HDF5 output for `Format=6`): build-time, not in the files.
