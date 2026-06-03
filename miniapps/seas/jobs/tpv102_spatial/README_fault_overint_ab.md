# TPV102 fault-flux over-integration A/B (Frontera)

Goal: measure **whether `--fault-overint` reduces the on-fault σ_n speckle/drift**
on the pure-upwind ADER **spatial** driver path (`seas_spatial_dyn_driver`).
Two runs, identical except one knob:

| run | flag | fault quadrature | expected on-fault σ_n |
|-----|------|------------------|------------------------|
| baseline | `--fault-overint 0` | degree `2·order` (=2 at p1) | speckles / drifts (TPV102: 120 → 80, −30, 75 MPa) |
| over-int | `--fault-overint 2` | degree `2·(order+2)` (=6 at p1) | should stay near the constant 120 MPa |

Both: **p1 · ADER-O2 · pure upwind · scalar Riemann**, same TOML
(`tpv102/configs/tpv102_spatial.toml`), same mesh, same `--deriv-cache
--shared-ck-recursion`. The only difference is the over-int factor — so any σ_n
change is attributable to over-integration alone.

Files:
- `tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch` — one parameterized job
  (factor via `FAULT_OVERINT`; output dir auto-encodes it so OFF/ON never collide).
- `submit_fault_overint_ab.sh` — chains build → baseline → over-int.

## 0. Prerequisite — be on the branch with the code

The `--fault-overint` knob is Phase 1+4 on `feat/fault-overint-resample`:
```
cd <seas-mfem repo>
git fetch origin
git checkout feat/fault-overint-resample
git log --oneline -2     # expect: 3521195 Phase 4 …  /  1cc0bbf Phase 1 …
```

## 1. One-shot submit (recommended)

```
# cheap first pass — 1000 m mesh, flex, 2 nodes, tfinal=5 (fast "does it help?")
CHEAP=1 bash miniapps/seas/jobs/tpv102_spatial/submit_fault_overint_ab.sh

# full production A/B — 200 m mesh, normal, 10 nodes/500 ranks, tfinal=12
bash miniapps/seas/jobs/tpv102_spatial/submit_fault_overint_ab.sh
```
The script submits a clean-rebuild **build** job, then the **baseline** and
**over-int** runs each with `--dependency=afterok` on the build. Override knobs:
`TPV102_TFINAL=<s>`, `NO_BUILD=1` (binary already built).

## 2. By hand (if you prefer)

```
# build once (the run jobs do NOT build):
sbatch miniapps/seas/jobs/spatial_dyn_build.sbatch clean-rebuild

cd miniapps/seas
FAULT_OVERINT=0 sbatch jobs/tpv102_spatial/tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
FAULT_OVERINT=2 sbatch jobs/tpv102_spatial/tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
```
Cheap 1000 m / flex variant of the same job:
```
FAULT_OVERINT=0 TPV102_MESH=tpv102/mesh/tpv102_1000m.msh TPV102_TFINAL=5 \
    sbatch -p flex -N 2 -n 100 -t 02:00:00 \
    jobs/tpv102_spatial/tpv102_p1_aderO2_pureupwind_fault_overint_ab.sbatch
# (and again with FAULT_OVERINT=2)
```

## 3. Confirm the knob actually fired

The over-int run prints, on rank 0:
```
[fault] over-integration ON (--fault-overint 2): fault-face quad degree 6 (baseline 2); QPs/face 12
```
Check it (absent ⇒ the run was the K=0 baseline):
```
grep -h "\[fault\] over-integration ON" tpv102_overint_ab_<ON_jobid>.out
```

## 4. Compare

Output dirs auto-encode the factor:
```
tpv102/out_p1_aderO2_pu_overint0_<jobid>/tpv102_station_*.dat   # baseline
tpv102/out_p1_aderO2_pu_overint2_<jobid>/tpv102_station_*.dat   # over-int
```
Plot/diff the **σ_n column** at each on-fault SCEC station. The plan's
prediction (§7.2): baseline drifts off 120 MPa; over-int stays much closer to
the analytically-constant 120 MPa, with `max_F |σ_n − σ_n0|` reduced and the
high-wavenumber speckle gone. Also sanity-check `V_max` is **bounded** in both.

## Notes / caveats

- **Pure upwind only.** `--fault-overint K>0` aborts if combined with mixed flux
  or precomputed face fluxes (Phase-1 scope). This job is pure upwind, so fine.
- **Byte-exact baseline.** `FAULT_OVERINT=0` is byte-identical to the pre-Phase-4
  driver (the driver skips `SetFaultOverint` at K=0), so the baseline IS the
  unmodified production result.
- **np>1 RS gate.** The driver flags parallel rate-state results as PRELIMINARY
  (shared-seam QPs use end-of-step ψ). Both runs see the same gate, so the A/B
  *comparison* is still valid; for a fully clean single-frame result run the
  cheap pass at smaller `-N`/`-n` or serially.
- **Matrix path (TPV31).** The same knob works on `interior_flux="matrix"`
  (the fault-flux routines live in base `WaveOperator`); a TPV31 A/B is the
  analogous job in `jobs/tpv31_spatial/` once you want the bi-material target.
