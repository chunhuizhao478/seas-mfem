# TPV102 fault-dealiasing A/B/C (Frontera)

Goal: measure **whether the two SeisSol dealiasing knobs reduce the on-fault σ_n
speckle/drift** on the pure-upwind ADER **spatial** driver path
(`seas_spatial_dyn_driver`).  Three runs, identical except the knobs:

| arm | flags | what it adds | expected on-fault σ_n |
|-----|-------|--------------|------------------------|
| A baseline | `--fault-overint 0` | minimal `2·order` fault quadrature (=2 at p1) | speckles / drifts (TPV102: 120 → 80, −30, 75 MPa) |
| B over-int | `--fault-overint 2` | degree `2·(order+2)`=6 (12 QPs/face): closes the **per-step** flux-aliasing channel | speckle gone; should stay near 120 MPa |
| C +resample | `--fault-overint 2 --fault-resample` | per-face degree-1 L2 projection of the per-macro-step **Δψ**: closes the **secular** state-accumulation channel | residual creep gone; flattest σ_n |

All three: **p1 · ADER-O2 · pure upwind · scalar Riemann**, same TOML
(`tpv102/configs/tpv102_spatial.toml`), same mesh, same `--deriv-cache
--shared-ck-recursion`.  The only differences are the two dealiasing knobs, so
any σ_n change is attributable to them.  Per plan §3.3/§4.2 the two are
**coupled**: resample is a no-op without over-integration (R = I unless
#QP > #DOF), so arm C requires `--fault-overint > 0` (the driver disables
resample otherwise).

Files:
- `tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch` — one parameterized job;
  knobs via `FAULT_OVERINT` and `FAULT_RESAMPLE` (output dir auto-encodes both so
  the three arms never collide).
- Standalone, self-contained variants (hardcoded knobs, no env, no dependency):
  `ab_overint0_*`, `ab_overint2_*`, `ab_overint2_resample_*` (each in a
  `_200m_normal` production and a `_1000m_flex` cheap form).
- `submit_fault_overint_ab.sh` — chains build → A → B → C.
- `compare_sigma_n.py` — 2- or 3-way σ_n-drift table (+ `--plot`).

## 0. Prerequisite — be on the branch with the code

The knobs are Phase 1-4 on `feat/fault-overint-resample`:
```
cd <seas-mfem repo>
git fetch origin
git checkout feat/fault-overint-resample
git log --oneline -3     # expect: …Phase 3 Δψ resample / Phase 2 R / …
```

## 1. One-shot submit (recommended)

```
# cheap first pass — 1000 m mesh, flex, 2 nodes, tfinal=5 (fast "does it help?")
CHEAP=1 bash miniapps/seas/jobs/tpv102_spatial/submit_fault_overint_ab.sh

# full production A/B/C — 200 m mesh, normal, 10 nodes/500 ranks, tfinal=12
bash miniapps/seas/jobs/tpv102_spatial/submit_fault_overint_ab.sh
```
The script submits a clean-rebuild **build** job, then the **baseline (A)**,
**over-int (B)** and **over-int+resample (C)** runs each with
`--dependency=afterok` on the build.  Override knobs: `TPV102_TFINAL=<s>`,
`NO_BUILD=1` (binary already built), `AB_NO_RESAMPLE=1` (original 2-arm only).

## 2. By hand (if you prefer)

Either the **standalone** files (hardcoded knobs, just `sbatch` them):
```
# build once (the run jobs do NOT build):
sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild
cd miniapps/seas
sbatch jobs/tpv102_spatial/ab_overint0_200m_normal.sbatch            # A
sbatch jobs/tpv102_spatial/ab_overint2_200m_normal.sbatch            # B
sbatch jobs/tpv102_spatial/ab_overint2_resample_200m_normal.sbatch   # C
# cheap 1000 m / flex forms: ab_*_1000m_flex.sbatch
```
or the **parameterized** job via env knobs:
```
cd miniapps/seas
FAULT_OVERINT=0                  sbatch jobs/tpv102_spatial/tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
FAULT_OVERINT=2                  sbatch jobs/tpv102_spatial/tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
FAULT_OVERINT=2 FAULT_RESAMPLE=1 sbatch jobs/tpv102_spatial/tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
# cheap 1000 m / flex: add TPV102_MESH=tpv102/mesh/tpv102_1000m.msh TPV102_TFINAL=5 -p flex -N 2 -n 100 -t 02:00:00
```

## 3. Confirm the knobs actually fired

Arm B prints, on rank 0:
```
[fault] over-integration ON (--fault-overint 2): fault-face quad degree 6 (baseline 2); QPs/face 12
```
Arm C prints that line **and**:
```
[fault] resample ON (--fault-resample): degree-1 L2 projector, 12 QPs/face
```
Check them (absent ⇒ that knob was off):
```
grep -h "\[fault\] over-integration ON" <B_jobid>.out
grep -h "\[fault\] resample ON"          <C_jobid>.out
```

## 4. Compare

Output dirs auto-encode the knobs:
```
tpv102/out_p1_aderO2_pu_overint0_<jobid>/tpv102_station_*.dat            # A baseline
tpv102/out_p1_aderO2_pu_overint2_<jobid>/tpv102_station_*.dat            # B over-int
tpv102/out_p1_aderO2_pu_overint2_resample_<jobid>/tpv102_station_*.dat   # C +resample
```
3-way σ_n-drift table (and PNGs):
```
cd miniapps/seas
python3 jobs/tpv102_spatial/compare_sigma_n.py <A_dir> <B_dir> <C_dir> --plot
```
The plan's prediction (§7.1/§7.2): A drifts off 120 MPa; **B** removes the
per-step speckle (`max_F |σ_n − σ_n0|` sharply down, checkerboard gone); **C**
additionally suppresses the residual secular creep (flattest σ_n, the `C/A`
reduction ≥ `B/A`).  Also sanity-check `V_max` is **bounded** in all three.

## Notes / caveats

- **Pure upwind only.** `--fault-overint K>0` aborts if combined with mixed flux
  or precomputed face fluxes (Phase-1 scope). This job is pure upwind, so fine.
- **Byte-exact baseline.** `FAULT_OVERINT=0` is byte-identical to the pre-Phase-4
  driver (the driver skips `SetFaultOverint` at K=0), so the baseline IS the
  unmodified production result.
- **Resample needs over-integration + rate-state.** `--fault-resample` is the
  *secular* layer; with over-int off it is disabled (a true no-op at every order,
  R-002).  TPV102 is rate-state (aging), so arm C resamples the Δψ increment.
  (LSW/TPV31 resamples the slip magnitude instead — same flag, different target.)
  Arm C's σ_n change vs B is attributable to resample alone (over-int is on in both).
- **np>1 RS gate.** The driver flags parallel rate-state results as PRELIMINARY
  (shared-seam QPs use end-of-step ψ). Both runs see the same gate, so the A/B
  *comparison* is still valid; for a fully clean single-frame result run the
  cheap pass at smaller `-N`/`-n` or serially.
- **Matrix path (TPV31).** The same knob works on `interior_flux="matrix"`
  (the fault-flux routines live in base `WaveOperator`); a TPV31 A/B is the
  analogous job in `jobs/tpv31_spatial/` once you want the bi-material target.
