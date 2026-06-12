# SAFS → SeisSol port — run guide

Files in this directory (ported from the MFEM SAFS LSW dynamic-rupture problem,
formatted to the updated **tpv13 / QuakeWorx** convention):

| File | Purpose |
|---|---|
| `parameters.par` | SeisSol namelist (FL=16 LSW, nucleation timing, output, EndTime) |
| `safs_fault.yaml` | easi fault model: friction, deep barrier, `Tnuc_s` patch; `!Include`s the stress file |
| `safs_initial_stress.yaml` | initial effective stress tensor (separated, tpv13 format) |
| `safs_material.yaml` | constant elastic medium (ρ/μ/λ) |

The physics: real San-Andreas segment, linear slip weakening (μ_s=0.85, μ_d=0.30,
D_c=2.0 m), regional H&Z effective background stress, time-ramped Gaussian strike
over-stress nucleation (20 MPa, 4 km, 1 s), deep barrier z∈[−20,−15] km.

**Format follows tpv13**: stress split into a `!Include`d file; `MeshFile` with the full
`.puml.h5` name; 11-field fault `OutputMask`; `numflux`+`numfluxnearfault`; `&Output` using
`Format`/`iPlasticityMask`/`Checkpoint`; empty `&Pickpoint`; `forced_rupture_time = 1e10` (disabled).

---

## Run

```bash
cd /path/to/this/dir
mkdir -p output
mpirun -np <N> /path/to/SeisSol_<config>  parameters.par
```
`OutputFile = 'output/safs'` requires the `output/` directory to exist first.

---

## Prerequisite: the mesh (NOT auto-generated here)

`parameters.par` points at `MeshFile = 'safs_mesh.puml.h5'`. Produce it from the SAFS gmsh
mesh (`safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh`):

1. **Retag boundaries to SeisSol codes** (`docs/fault-tagging.rst`): MFEM uses
   fault=101, top=102, bottom=103, sides=104, rock=1. SeisSol needs:

   | MFEM tag | SeisSol boundary code |
   |---|---|
   | 101 fault | **3** (dynamic rupture) |
   | 102 top | **1** (free surface) |
   | 103 bottom | **5** (absorbing) |
   | 104 sides | **5** (absorbing) |
   | 1 rock | volume group (kept) |

2. **Convert with pumgen**: `pumgen safs_mesh.msh -s msh2` → `safs_mesh.puml.h5`.

   ⚠️ `pumgen` is not installed locally (prior session) — build it or run the conversion on a
   cluster. The param files do not depend on having pumgen locally.

3. SeisSol requires a **right-handed** coordinate system; the fault must not lie in the xy-plane
   (`docs/dynamic-rupture.rst`) — the near-vertical SAF satisfies this.

---

## Build notes (not expressible in the .par)

- **Order**: MFEM `ader_order = 2` maps to SeisSol's compile-time `CONVERGENCE_ORDER` (build with
  order 2 to match MFEM, or ≥3 for fidelity).
- **Flux**: MFEM `mixed_flux = none` (pure upwind) == SeisSol `numflux='godunov'`.
- **HDF5**: `Format=6` (hdf5 wavefield) needs a SeisSol HDF5 build; set `Format=10` to disable
  volume output (recommended first on the large SAFS mesh).

---

## Validation gates (first SeisSol run)

### Initial-traction signs / stress port
From the fault output read `Ts0`,`Td0`,`Pn0` at the hypocenter ≈ (606971, 3707270, −4965.62).
MFEM reports (comp-positive): σ_n_eff ≈ 49.27 MPa, τ_strike0 ≈ 26.32, τ_dip0 ≈ −12.65, |τ_pre| ≈ 29.2.
- `Pn0` should be compressive (< 0) with |Pn0| ≈ 49 MPa.
- If `|Pn0|` matches but `Ts0`/`Td0` flip sign across the fault → the strike/dip frame is mirrored;
  **fix via `XRef/YRef/ZRef`** (REVIEW.md R-001), not by editing stresses.

### Nucleation
Rupture should initiate inside the 4 km patch during the first ~1 s (the `Tnuc_s` smoothStep ramp),
then propagate spontaneously. If it drives the wrong slip sense, flip the sign of `Tnuc_s`.

---

## Intentional differences from MFEM (do not "fix")

- **No σ_n strength floor / "normalcap" (10 MPa)**: an MFEM numerical stabilizer with no SeisSol
  equivalent (SeisSol LSW uses `min(σ_n,0)`). If instability appears, SeisSol's lever is `etaDamp`
  (`&DynamicRupture`, <1.0).
- **Pore pressure** is baked into the effective stress diagonal (SeisSol LSW does not subtract P_p).
- **`forced_rupture_time = 1e10`** (disabled): nucleation is purely the `Tnuc_s` over-stress.
  Including it (vs omitting) matches the tpv13 style and is safe — with 1e10 the forced-rupture
  term stays 0 and does not affect the `Tnuc_s`/`s_0`/`t_0` ramp.
- **No off-fault plasticity** (`Plasticity = 0`): MFEM SAFS had none, so `safs_material.yaml` omits
  tpv13's `plastCo`/`bulkFriction` and the volume stress `!Include`.
- **No fault receivers**: `OutputPointType = 4` (ParaView only). To add ASCII receivers, set
  `OutputPointType = 5`, add `PPFileName` to `&Pickpoint`, and a receiver coordinate file.
