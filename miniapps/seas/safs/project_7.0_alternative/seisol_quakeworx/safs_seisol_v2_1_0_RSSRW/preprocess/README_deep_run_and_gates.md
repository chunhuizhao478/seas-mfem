# RSSRW on `safv4_deep_500m` — run guide + revised validation gates

Companion to `PLAN_safv4_deep_mesh_port_2026-06-15.md` (Phases 3–5). This is the
RSSRW physics (FL=103 rate-and-state + strong velocity weakening, CVM material via
ASAGI, Gaussian over-stress nucleation) run on the 3-strand `safv4_deep_500m` mesh
instead of the single-fault baseline. The physics inputs are reused VERBATIM; only
the mesh changed.

## Status (as built)

- Mesh retagged + converted + verified (Phases 0–2, PASS). The single DR group is
  faults 101/102/103 → BC 3; top 201 → BC 1; bottom 202 + sides 203 → BC 5.
  PUML check: tets 1,180,159 / nodes 279,448; BC1 175,033; BC5 122,593; BC3 264,666
  (= 2× the 132,333 fault triangles — this pumgen tags both adjacent tet faces).
- `parameters.par` now points at `MeshFile = 'safs_mesh_deep.puml.h5'` (in place).
- Reused unchanged: `safs_fault.yaml`, `safs_initial_stress.yaml`,
  `safs_material_cvm.yaml`, `safs_material_cvm.nc`.

## Run (cluster or QuakeWorx)

Final-submission file set (all at the TOP level of `safs_seisol_v2_1_0_RSSRW/`,
NOT `preprocess/`):
```
parameters.par
safs_fault.yaml
safs_initial_stress.yaml
safs_material_cvm.yaml
safs_material_cvm.nc            (112 MB; use the gateway data manager if upload-limited)
safs_mesh_deep.puml.h5         (the deep mesh)
```
The single-fault `safs_mesh.puml.h5` is NOT part of this run — do not upload it.

```bash
cd /path/to/safs_seisol_v2_1_0_RSSRW
mkdir -p output            # OutputFile='output/safs' requires output/ to exist first
mpirun -np <N> /path/to/SeisSol_<order>_<...>  parameters.par
```
Build: CONVERGENCE_ORDER ≥ 2 (MFEM ader_order=2), ASAGI enabled (the QuakeWorx app
has it). `numflux='godunov'` = MFEM pure upwind. `Format=10` (no volume output) is set
for the large mesh; fault `&Elementwise`, free-surface, and energy outputs are on.

Smoke first: optionally drop `EndTime` (`&AbortCriteria`) to ~10–20 s to reach the
time loop and first propagation cheaply, then restore `EndTime = 200.0` for production.

---

## Validation gates

### Gate A — init reaches the time loop
No easi/ASAGI exception. Log must show the RS parameter source:
```
RS parameter source (1 == from parameter file, 0 == from easi file): f0 1 - muW 1 - b 0
```
`b` MUST be `0` (spatial `rs_b` from easi). A constant-`b` run is WRONG away from the
hypocenter (`safs_fault.yaml` header; the `RS_b=0.0168` in `parameters.par` is a
fallback only).

### Gate B — material / CVM coverage (REVISED for the deep mesh)
**ASAGI `out of range. Fixing.` warnings are now EXPECTED — the old "zero out-of-range"
gate no longer applies.** The deep mesh extends beyond the regional CVM grid
(`x[303000,697500] y[3612000,3903000] z[-45000,0]`):

- **39.2%** of mesh nodes fall outside the grid xy-footprint (far-field),
- **28.1%** of nodes are above z=0 (real topography, up to +3006 m),
- 0% below z=−45000.

ASAGI clamp-to-boundary supplies edge / z=0 values there. This is acceptable: the
far-field is the absorbing region, and clamping near-surface velocity is reasonable.

New gate wording: clamp warnings must occur ONLY outside the documented footprint /
above z=0. The near-fault region IS inside the grid and must be un-clamped — spot-check
the material at the hypocenter facet:
- **μ ≈ 23.4 GPa** (the CVM value; `safs_fault.yaml` / CALC doc). NOT the constant
  fallback μ=3.2e10 (`safs_material_cvm.yaml` `!ConstantMap`) and NOT a clamped edge
  value.

Do NOT regenerate `safs_material_cvm.nc`: the converter rejects `z_max`/xy beyond the
regional sidecar footprint (`convert_cvm_to_asagi.py` z-range guard), so the grid
cannot be extended to cover the topography/far-field. Reuse it as-is.

A/B fallback still available: set `MaterialFileName = 'safs_material.yaml'` (constant
medium) in `parameters.par` for debugging — that file lives in `preprocess/archive/`,
move it to the top level if you use it.

### Gate C — initial tractions / stress port
At the hypocenter facet read `Ts0`,`Td0`,`Pn0`:
- `Pn0` compressive (< 0), `|Pn0|` ≈ 49 MPa.
- If `|Pn0|` matches but `Ts0`/`Td0` flip sign across the fault → the strike/dip frame
  is mirrored; fix via `XRef/YRef/ZRef` (REVIEW R-001), NOT by editing stresses.

### Gate D — strike/dip frame across strands (NEW, multi-strand specific)
The single reference `XRef/YRef/ZRef = (0,-1,0)`, `refPointMethod=1` was validated for
the single near-vertical SAF strand. The deep mesh faults are **moderately dipping**
(Phase-0 finding: median |n_z| ≈ 0.44 ⇒ median dip ≈ 64°, ~12% of facets dip <46° in
the San Gorgonio Pass geometry) and there are THREE differently-oriented strands. So:
- Confirm each strand's `Ts0/Td0` signs are physical (right-lateral on the SAF).
- If a strand comes out mirrored, iterate on `XRef/YRef/ZRef` only and document the
  change here. This is the main NEW physics risk introduced by the geometry.

### Gate E — nucleation
The `Tnuc_s` Gaussian (24 MPa, R=6000 m, centre (606971, 3707270, −4965.62), smoothStep
over [0,1] s) seeds rupture in the first ~1 s, then it propagates spontaneously.
**Phase-0 finding:** the patch centre is closest to the **Garnet Hill strand
(12.6 m)**, not the main San Andreas — nucleation seeds in the Garnet Hill region, but
since all strands are one DR group the 6 km Gaussian also covers nearby SAF facets.
Confirm rupture initiates inside the patch and spreads onto the SAF.

### Gate F — stability / sanity
No NaN, no dt→0, no MFEM/SeisSol abort. (A `PMPI_Waitall "request pending due to
failure"` with low MaxRSS and no NaN on Frontera is a fabric blip, not our bug —
resume from checkpoint; but `Checkpoint=0` here, so a smoke just re-runs.)

---

## Provenance

- Source mesh: `project_7.0_preferred/meshing/results/safv4_deep_500m.msh`
  (1,180,159 tets; strands 101 San Andreas / 102 Banning / 103 Garnet Hill;
  boundaries 201 top / 202 bottom / 203 sides).
- Retag: `preprocess/retag_safv4_deep_to_seissol.py` → `preprocess/safs_mesh_deep.msh`.
- Convert: `preprocess/build_and_run_pumgen_linux.sh` (Linux) → `safs_mesh_deep.puml.h5`.
- Verify: `preprocess/verify_puml_deep.py` (PASS) and `preprocess/preflight_deep_mesh.py`
  (PASS).
