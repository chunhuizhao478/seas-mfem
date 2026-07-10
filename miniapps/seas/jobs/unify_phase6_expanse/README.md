# Unify-plan Phase 6 — Expanse verification + Caliper speed jobs

Production-mesh verification of the **unified interior/shared fault substep path**
(`PLAN_unify_interior_shared_fault_substep_2026-07-09.md`, Phase 6 steps 1+3) plus
per-polynomial-order speed tracking, on the branch
`system/substep-unified-tpv26-tpv27`.

## Jobs (six, one per problem × polynomial order)

| file | problem | law | P / ADER-O | wall |
|---|---|---|---|---|
| `tpv205_p1_npsweep_cali_expanse.sbatch` | TPV205 (200 m) | LSW | 1 / 2 | 12 h |
| `tpv205_p2_npsweep_cali_expanse.sbatch` | | | 2 / 3 | 24 h |
| `tpv205_p3_npsweep_cali_expanse.sbatch` | | | 3 / 4 | 48 h |
| `tpv102_p1_npsweep_cali_expanse.sbatch` | TPV102 (200 m) | RS aging | 1 / 2 | 12 h |
| `tpv102_p2_npsweep_cali_expanse.sbatch` | | | 2 / 3 | 24 h |
| `tpv102_p3_npsweep_cali_expanse.sbatch` | | | 3 / 4 | 48 h |

Each job runs three phases:

- **A (verify)** — np sweep (default `16 32 64 128`, one node) of
  `--fault-iterator substep`, short horizon (default 0.5 s), station traces only
  (no ParaView, no checkpoints).
- **B (compare)** — `compare_np_sweep.py` prints max |X_np − X_ref| per station
  (ref = smallest np in the sweep) into the log and
  `<outdir>/RESULT_np_sweep.txt`; non-zero exit if any station exceeds
  `UNIFY_CMP_TOL` (default 1e-9, the trace print precision).  With the unified
  path + the deterministic station tie-break
  (`np4_attractor_root_cause_2026-07-10.md`), the expectation is print-precision
  equality; a violation is the plan's Phase-6 STOP condition.
- **C (Caliper)** — one profiled run per job (default 128 ranks, 1.2 s) with
  `runtime-report` + `spot` channels, written as
  `{prob}_p{P}_cali_<jobid>.cali-region-report.txt` / `.cali` next to the job
  logs.  Skipped with a warning if the binary is not Caliper-linked.

## Build (once, on Expanse)

```bash
USE_CALIPER=YES bash <seas-mfem-root>/build_expanse.sh
cd <seas-mfem-root>/miniapps/seas
make seas_tpv205_driver seas_tpv102_driver -j 16
```

Meshes are gitignored — generate from the committed `.geo` (Gmsh **v2.2 ASCII**,
see CLAUDE.md's msh22 requirement) or scp them:

```bash
gmsh tpv205/mesh/tpv2053d_200m.geo -3 -format msh22 -o tpv205/mesh/tpv2053d_200m.msh
gmsh tpv102/mesh/tpv102_200m.geo   -3 -format msh22 -o tpv102/mesh/tpv102_200m.msh
```

## Submit

```bash
cd <seas-mfem-root>/miniapps/seas/jobs/unify_phase6_expanse
sbatch tpv205_p1_npsweep_cali_expanse.sbatch     # etc.
```

Knobs (via `sbatch --export=ALL,VAR=...`): `TPV_NP_SWEEP`, `TPV_MESH`,
`TPV_TFINAL_VERIFY`, `TPV_TFINAL_CALI`, `TPV_CALI_RANKS`, `TPV_CFL`,
`TPV_SKIP_VERIFY`, `TPV_SKIP_CALI`, `TPV_CALI_CONFIG`, `TPV_OUTDIR`,
`UNIFY_CMP_TOL`.  Note: `np ≤ 8` on the 200 m meshes will not fit the wall
time — only use a small-np sweep together with a coarser `TPV_MESH`
(e.g. `tpv102/mesh/tpv102_1000m.msh`, sweep `"1 2 4 8"`), which is also the
only way to anchor the sweep at np=1 on Expanse.  The np=1 anchor for the
200 m sweeps is established locally (symmirror + smoke evidence in
`PHASE6_sym1000_revalidation_2026-07-10.md`).

## After all six jobs finish

```bash
python3 plot_caliper_speed.py         # -> caliper_speed_vs_order.{csv,png}
```

produces the speed-vs-polynomial-order graph (avg time/rank per headline
region, max-rank caps showing load imbalance) for both problems, plus the CSV
for the plan's Phase-6 record.  Wall-time guesses in the `#SBATCH --time`
headers are deliberately generous — trim them after the first P1 timings.
