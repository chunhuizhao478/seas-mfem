# SeisSol TPV104 — paired with MFEM full-domain mesh

This directory ships a complete SeisSol TPV104 input set that produces
output station files directly comparable to MFEM's `seas_tpv104_driver`
output.  The setup is a port of the canonical
[SeisSol/Examples/tpv104](https://github.com/SeisSol/Examples/tree/master/tpv104)
(local mirror at `/Users/chunhuizhao/projects/seisol_exmaples/tpv104/`)
with **one structural deviation**: a **full-domain unstructured tet
mesh** instead of upstream's half-mesh-and-mirror.  Physics parameters
(RS_muW, t_0, RS_b, etc.) are byte-identical to upstream.

## Why not the upstream half-mesh-and-mirror?

Upstream uses `tpv104_half.geo` (only y ≥ 0), then
[mirrorMesh.py](https://github.com/SeisSol/Meshing/blob/master/mirrorMesh/mirrorMesh.py)
reflects across y = 0 to produce `tpv104_half_sym.h5`.  This guarantees
the resulting mesh is exactly mirror-symmetric across the fault plane,
eliminating one source of σ_n perturbation discretisation residue
(see `debug_document/tpv104_debug_document/tpv104_mesh_asymmetry_finding_2026-04-24.md`).

For our cross-verification with MFEM, both codes must run on the
**same** full-domain unstructured mesh so any MFEM-vs-SeisSol delta is
attributable to numerics, not to mesh-symmetry differences.  Therefore
this directory uses MFEM's `tpv104/mesh/tpv104_*m.geo` byte-for-byte,
modulo a relabelling of `Physical Surface` tags (1/3/5 → 101/103/105)
to match SeisSol's PUMGen convention.

## File inventory

| file | purpose |
|---|---|
| `tpv104_seisol_500m.geo` | Gmsh source, 500 m fault resolution; identical geometry to MFEM `tpv104_500m.geo` modulo Physical Surface tag remap |
| `tpv104_seisol_1000m.geo` | 1000 m coarse dev variant |
| `tpv104_seisol_material.yaml` | Homogeneous half-space (ρ, μ, λ); byte-identical to upstream |
| `tpv104_seisol_fault.yaml` | Initial stress + per-QP friction (rs_a, rs_srW, rs_sl0); LuaMap port of upstream |
| `tpv104_seisol.par` | Main parameters; FL=103, RS_muW=0.2, t_0=0.5 (upstream values) |
| `tpv104_seisol_faultreceivers.dat` | 9 on-fault stations matching SCEC TPV104 / MFEM |
| `tpv104_seisol_receivers.dat` | 6 free-surface stations matching SCEC TPV104 / MFEM |
| `tpv104_seisol_mesh_build.sbatch` | gmsh + pumgen → PUMGen .h5/.xdmf |
| `tpv104_seisol_t12_run.sbatch` | t = 12 s production run, 400 ranks |

## Critical finding: MFEM's `f_w` is wrong

Cross-referencing SCEC's `benchmark_document/SCEC_validation_slip_law.pdf` table
against `seisol_exmaples/tpv104/parameters.par`:

| parameter | SCEC benchmark | SeisSol upstream | MFEM `config/tpv104_params.hpp` | status |
|---|---|---|---|---|
| `f_w` (= `RS_muW`) | **0.2** | **0.2** | **0.1** | ❌ **MFEM bug — half the spec value** |
| `T_nuc` (= `t_0`) | **0.5 s** (per upstream) | **0.5 s** | **1.0 s** | ❌ MFEM doubles upstream's |
| `s_xy` sign | (Cartesian) | `s_xy = +40e6` | `tau2_0 = +40e6` (BP5 t1=dip / t2=strike) | ✓ equivalent under their respective frames |
| `nuc_xy` sign | (Cartesian) | `nuc_xy = +45e6·F(r)` | `tau2_nuc += +Δτ₀·F(r)` (BP5 frame) | ✓ same physical perturbation |

**Implication**: MFEM's V_strike over-shoot (1.7-2.7×) vs the SeisSol reference
traces in `benchmark_data/seisol/` is **partly explained by an MFEM bug**:
`config/tpv104_params.hpp` sets `f_w = 0.1`, half the SCEC benchmark value of
0.2.  Lower `f_w` → friction collapses MORE at the rupture front → larger V.
Same direction as the over-shoot we observed.

Recommended fix in `config/tpv104_params.hpp`:

```cpp
-static constexpr real_t f_w    = 0.1;            ///< Weakening friction coefficient
+static constexpr real_t f_w    = 0.2;            ///< Weakening friction coefficient (SCEC TPV104 benchmark)
-static constexpr real_t nuc_T  = 1.0;            ///< Rise time [s]
+static constexpr real_t nuc_T  = 0.5;            ///< Rise time [s] (SCEC TPV104 + SeisSol Examples upstream)
```

The `f_w = 0.1` in MFEM is a misread of the SCEC PDF (the value that appears
elsewhere in the TPV104 documentation is the **friction coefficient at static
locking**, not the velocity-weakening asymptote `f_w`).  Update the
`Why:` comment chain in `tpv104_review_*` docs accordingly.

After the fix lands, rerun the t12 production sbatch and re-evaluate
the V_strike over-shoot.  If it drops substantially, parameter-misalignment
WAS the dominant cause of the over-shoot, and the wave-operator /
mesh-asymmetry investigation can be deprioritised.

Track in
`debug_document/tpv104_debug_document/tpv104_v_overshoot_divergence_2026-04-24.md`.

## Submission flow on Frontera

```bash
# 1. Copy the seisol/ directory to a Frontera scratch location.
cd /scratch2/$USER/tpv104_seisol_run
cp -r ${SEAS_MFEM}/miniapps/seas/jobs/tpv104/seisol/* .

# 2. Build the mesh (15 min on dev queue).
sbatch tpv104_seisol_mesh_build.sbatch
# Wait for completion, verify tpv104_seisol_500m.h5 + .xdmf are present.

# 3. Submit the production run.  Estimated 6-10 h wall-clock.
sbatch tpv104_seisol_t12_run.sbatch

# 4. After completion, copy `output/tpv104_seisol-faultreceiver-*` files
#    back to local and run the visualizer:
#       python3 ${SEAS_MFEM}/miniapps/seas/tpv104/visualize_results.py \
#           --mfem ${MFEM_RESULTS_DIR} --seisol --save
#    NOTE: the visualizer currently expects SeisSol output in
#    `benchmark_data/seisol/` flat-file layout.  An adapter is needed
#    to ingest fresh SeisSol-format `output/tpv104_seisol-faultreceiver-*`
#    runs.  TBD next session.
```

## Mesh-symmetry observation (relevant for σ_n leak diagnostic)

Even with the full-domain unstructured mesh, gmsh's Delaunay
tetrahedralisation produces a non-mirror-symmetric mesh — the same one
MFEM uses.  See
`debug_document/tpv104_debug_document/tpv104_mesh_asymmetry_finding_2026-04-24.md`
for the verification.  This means **both** SeisSol and MFEM, on this
full-domain mesh, will exhibit some level of σ_n discretisation
residue (≤ a few MPa at peak).  The σ_n perturbation under SeisSol on
this mesh is the **upper bound** on what MFEM should display once the
muW + t_0 parameters are aligned; differences above that bound are
attributable to other factors.

## Pre-flight checklist before submission

- [ ] User approves Frontera sbatch (per `feedback_frontera_approval.md`).
- [ ] SeisSol Frontera build verified at the path in
  `tpv104_seisol_t12_run.sbatch:SEISSOL_BIN`.  Build it per the
  `reference_seissol_frontera.md` memory (ICX 2023.1 + IMPI 2021.9 +
  ORDER=5) if not present.
- [ ] PUMGen built and on `$PATH` on Frontera.
- [ ] User decides muW/t_0 — match SeisSol upstream (this directory's
  defaults) OR match SCEC PDF (re-edit `tpv104_seisol.par` + flag in
  the debug document).

## Troubleshooting

* **PUMGen "missing boundary tag"**: verify the .geo's `Physical Surface`
  block uses 101/103/105.  The header sed transformation in
  `tpv104_seisol_500m.geo` should have set this; check by
  `grep "Physical Surface" tpv104_seisol_500m.geo`.
* **SeisSol parameters.par parse error**: SeisSol's Fortran namelist
  parser is finicky about `&block / /` syntax — preserve the trailing
  `/` on each block.
* **"easi: rs_a not found"**: confirm `tpv104_seisol_fault.yaml` is in
  the run cwd (not just the seisol/ source dir).
