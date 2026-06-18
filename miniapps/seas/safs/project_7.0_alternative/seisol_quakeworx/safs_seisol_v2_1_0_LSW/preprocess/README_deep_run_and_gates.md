# LSW on `safv4_deep_500m` — run guide + revised validation gates

Companion to `PLAN_safv4_deep_mesh_port_2026-06-15.md` (Phases 3–5; see the "LSW twin
adaptation note" box at the top of that file). This is the LSW physics (FL=16 linear
slip-weakening, CVM material via ASAGI, Gaussian over-stress nucleation) run on the
3-strand `safv4_deep_500m` mesh instead of the single-fault baseline. The physics
inputs are reused VERBATIM from `safs_seisol_v2_0_0_LSW`; only the mesh changed
(`safs_mesh.puml.h5` → `safs_mesh_deep.puml.h5`).

This is the LSW twin of `safs_seisol_v2_1_0_RSSRW` — same deep mesh, slip-weakening
friction instead of rate-and-state + strong velocity weakening.

## Status (as built)

- Mesh retagged + converted + verified (Phases 0–2, PASS). The single DR group is
  faults 101/102/103 → BC 3; top 201 → BC 1; bottom 202 + sides 203 → BC 5.
  PUML check (`verify_puml_deep.py`, re-run in this folder): tets 1,180,159 /
  nodes 279,448; BC1 175,033; BC5 122,593; BC3 264,666 (= 2× the 132,333 fault
  triangles — this pumgen tags both adjacent tet faces).
- The deep mesh artifacts (`safs_mesh_deep.msh`, `safs_mesh_deep.puml.h5`) are
  byte-identical to the RSSRW deep port; the verified PUML was COPIED here rather
  than re-running pumgen (the mesh is friction-agnostic).
- `parameters.par` points at `MeshFile = 'safs_mesh_deep.puml.h5'` (in place); only
  the `&MeshNml` block changed vs `safs_seisol_v2_0_0_LSW` — the friction
  (`FL = 16`), `OutputMask`, nucleation, and `EndTime` are LSW, unchanged.
- Reused unchanged from the LSW baseline: `safs_fault.yaml` (LSW model),
  `safs_initial_stress.yaml`, `safs_material_cvm.yaml`, `safs_material_cvm.nc`.

## Run (cluster or QuakeWorx)

Final-submission file set (all at the TOP level of `safs_seisol_v2_1_0_LSW/`,
NOT `preprocess/`):
```
parameters.par
safs_fault.yaml                (LSW: mu_d/d_c/cohesion + mu_s deep barrier + Tnuc_s)
safs_initial_stress.yaml
safs_material_cvm.yaml
safs_material_cvm.nc            (112 MB; use the gateway data manager if upload-limited)
safs_mesh_deep.puml.h5         (the deep mesh)
```
The single-fault `safs_mesh.puml.h5` is NOT part of this run — do not upload it.

```bash
cd /path/to/safs_seisol_v2_1_0_LSW
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

### Gate A — init reaches the time loop (LSW form)
No easi/ASAGI exception. SeisSol reports `Friction law 16` (linear slip-weakening) and
reads the LSW model from easi:
- `[mu_d, d_c, cohesion]` (ConstantMap: 0.12 / 1.8 m / 0.0),
- `[mu_s]` (LuaMap: the deep barrier — `1.0e6` for `z ∈ (−20000, −15000)`, `0.43`
  elsewhere),
- `[Tnuc_s]` (LuaMap Gaussian) and `[Tnuc_n, Tnuc_d, forced_rupture_time]`.

There is NO "RS parameter source" log line here — that is an FL=103 (rate-and-state)
diagnostic and does not apply to LSW. (The RSSRW twin checks `b 0`; the LSW twin has
no analogous check.)

### Gate A′ — deep barrier vs the deep mesh (NEW, LSW-specific; verified locally)
The LSW deep barrier is a fixed `z`-band, so it must be checked against the deep mesh's
fault depth extent. `preflight_deep_mesh.py` on the source mesh + a direct parse of the
retagged fault group give:
- fault z-extent **−19,329 … +2,210 m** (deepest facet **−19.3 km**),
- **0 fault facets below z = −20,000 m**, 17.3 % below −15,000 m.

⇒ the barrier `z ∈ (−20000, −15000)` correctly locks the deepest ~4.3 km of fault
(−15 → −19.3 km) and no facet falls below the band's lower edge, so nothing is left
ruptable below the barrier. **No barrier edit is needed** — this is why the LSW
`safs_fault.yaml` is reused verbatim. (At run time, confirm the deepest fault facets
show `mu_s = 1e6` / remain locked and that slip stays above ~15 km depth in the
seismogenic zone.)

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
At the hypocenter facet read `Ts0`,`Td0`,`Pn0` (same `safs_initial_stress.yaml` as the
baseline, so the same expectation):
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
  change here. This is the main NEW physics risk introduced by the geometry. (Identical
  to the RSSRW twin's Gate D — it is mesh/geometry-driven, not friction-driven.)

### Gate E — nucleation (LSW form)
The `Tnuc_s` Gaussian (**30 MPa** — 2026-06-15: raised 10 → 30 MPa, ~RSSRW's 24 MPa;
R=6000 m, centre (606971, 3707270, −4965.62), smoothStep over [0,1] s) seeds rupture in
the first ~1 s, then should propagate spontaneously. `Tnuc_n = Tnuc_d = 0` (strike-only)
and `forced_rupture_time = 1e10` (forced rupture disabled; nucleation is purely the
`Tnuc_s` ramp).

**Phase-0 finding:** the patch centre is closest to the **Garnet Hill strand
(12.6 m)**, not the main San Andreas — nucleation seeds in the Garnet Hill region, but
since all strands are one DR group the 6 km Gaussian also covers nearby SAF facets.
Confirm rupture initiates inside the patch and spreads onto the SAF.

**Why 30 MPa (deep-mesh smoke run, 2026-06-15):** the original **10 MPa** ARRESTED
immediately on the deep mesh — Mw 3.8, seismic moment plateaued at ~6.3e14 N·m. The
diagnostic was conclusive: peak fault slip was only **0.017 m = 0.9 % of Dc=1.8 m**
(0 of 264,666 cells reached Dc), so the slip-weakening law never activated → no stress
drop → no propagation. RSSRW propagates from the SAME Garnet Hill seed at 24 MPa, so the
fix is amplitude, not geometry. 30 MPa lifts peak shear at the seed facet to
≈51 MPa (>> strength 27.7 MPa). **This is the new test** — confirm whether the patch now
slips past Dc (≈1.8 m) and launches a front. The earlier `safs_fault.yaml` "candidate
D" arrest history (B at 22.8 km, C at +42 km in the s = 40–70 km low-`mu_apparent`
corridor) was on the single-fault mesh and is the *secondary* arrest risk if the front
does launch.

If 30 MPa still arrests, next levers (in order): `forced_rupture_time` in the patch
(standard SCEC LSW nucleation — guarantees the patch weakens to `mu_d`), then `d_c 1.8 →
1.4` + a finer fault mesh (Dc=1.8 m is forced by the 500 m cohesive-zone resolution,
Lambda_0 = 1890 m = 3.8 elements). NOT a mesh boundary/coordinate edit — the geometry is
confirmed good.

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
- LSW twin: the deep mesh (`.msh` + `.puml.h5`) and the deep-mesh tooling are reused
  byte-for-byte from `safs_seisol_v2_1_0_RSSRW`; the only LSW-specific inputs are
  `parameters.par` (`FL = 16`, mesh line) and `safs_fault.yaml` (LSW model), both
  carried over from `safs_seisol_v2_0_0_LSW`. The RSSRW-only friction calc
  (`CALC_rssw_friction_lengthscales`) is intentionally not reproduced; the LSW friction
  history lives inline in `safs_fault.yaml`.
