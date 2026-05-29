# Caliper "perfgraph" — build & run runbook (TPV31 spatial dynamic-rupture)

This is the operator-facing recipe for producing a per-region timing tree
("perfgraph") from `seas_spatial_dyn_driver`. It accompanies
`caliper_perfgraph_plan.md` (design) / `_check.md` (review) / `_fix.md`.

## What was instrumented

`MFEM_PERF_SCOPE` / `MFEM_PERF_BEGIN`/`MFEM_PERF_END` regions
(`general/annotation.hpp`) were added to the dynamic-rupture hot path. On a
**default build (`MFEM_USE_CALIPER=NO`) these are no-ops** — zero numerical and
zero performance impact. They only emit data on a Caliper-enabled build.

Regions emitted (the perfgraph node names):

| Region | Path | Location |
|--------|------|----------|
| `seas::spatial_dyn::step` | both | driver per-macro-step loop body |
| `seas::spatial_dyn::AdvanceADERWithSubStep` | ADER | driver macro-step orchestration |
| `seas::spatial_dyn::friction_substep` | both | per-QP LSW/RS friction solve (`iterator.Advance`) |
| `seas::WaveOperator::ComputeADERSubStepStates` | ADER | ADER Cauchy–Kovalewski predictor |
| `seas::WaveOperator::ApplySpatialDerivative` | ADER | predictor per-direction spatial derivative |
| `seas::WaveOperator::AdvanceADER` | ADER | ADER corrector |
| `seas::WaveOperator::ComputeADERVolumeUpdate` | ADER | ADER volume update (→ `ComputeVolumeRHS`) |
| `seas::WaveOperator::ComputeADERFaceFluxRHS` | ADER | ADER interior-face flux assembly |
| `seas::WaveOperator::ComputeADERSharedFaceFluxRHS` | ADER | ADER MPI-seam face flux assembly |
| `seas::WaveOperator::ComputeVolumeRHS` | both | volume term (shared by ADER + RK) |
| `seas::WaveOperator::Mult` | RK only | semidiscrete RHS (`--time-integrator rk4\|rk45`) |
| `seas::WaveOperator::ComputeFaceFluxRHS` | RK only | interior-face flux (via `Mult`) |
| `seas::WaveOperator::ComputeSharedFaceFluxRHS` | RK only | MPI-seam face flux (via `Mult`) |
| `seas::WaveOperator::ComputeMaxDt` / `seas::BimaterialWaveOperator::ComputeMaxDt` | both | CFL dt estimate (once, setup) |

The expected nesting for the **ADER path** (TPV31; `--time-integrator ader`):

```
step
└─ AdvanceADERWithSubStep
   ├─ ComputeADERSubStepStates   (→ ApplySpatialDerivative)
   ├─ friction_substep
   └─ AdvanceADER
      ├─ ComputeADERVolumeUpdate (→ ComputeVolumeRHS)
      ├─ ComputeADERFaceFluxRHS
      └─ ComputeADERSharedFaceFluxRHS
ComputeMaxDt                      (once, at setup — top level)
```

So friction-solve time (`friction_substep`) is readable separately from the ADER
predictor/corrector and its volume vs. face-flux assembly. The `Mult` /
`ComputeFaceFluxRHS` / `ComputeSharedFaceFluxRHS` regions fire ONLY on the
explicit-RK path (`--time-integrator rk4|rk45`) and are **absent from a TPV31
(ADER) perfgraph** — TPV31's face-flux cost shows under
`ComputeADERFaceFluxRHS` / `ComputeADERSharedFaceFluxRHS` instead.

> **Per-face hooks (`InteriorFaceFlux_` / `SharedInteriorFaceFlux_` /
> `ApplyElementJacobian_`) are deliberately NOT annotated** — they are called
> O(faces) times per step; their cost is captured by the enclosing
> `Compute{,ADER}FaceFluxRHS` / `Compute{,ADER}SharedFaceFluxRHS` regions
> without per-call Caliper overhead.

## Prerequisites (BOTH required; neither is the default)

1. **MFEM + the seas driver built with `MFEM_USE_CALIPER=YES`** (links
   `libcaliper`). The committed default is `MFEM_USE_CALIPER = NO`
   (`config/defaults.mk:181`) — do not commit `=YES`.
2. **These annotations present** (this branch). On any build without (1) the
   `TPV31_PERFGRAPH=1` sbatch hook is inert and prints nothing useful.

## Build recipe (Frontera)

```bash
# 1. Obtain Caliper (once). Record the install prefix as CALIPER_DIR.
spack install caliper +adiak          # or build from github.com/LLNL/Caliper (CMake)
CALIPER_DIR=$(spack location -i caliper)

# 2. Reconfigure + rebuild MFEM with Caliper (in the MFEM build dir).
make config MFEM_USE_CALIPER=YES CALIPER_DIR="$CALIPER_DIR"   # [ADIAK_DIR=... optional]
make -j

# 3. Rebuild the seas driver — FORCE-CLEAN.  The seas Makefile has no header
#    deps, so an incremental rebuild after flipping Caliper would relink STALE
#    objects that compiled every MFEM_PERF_* as a no-op → libcaliper links but
#    ZERO seas regions appear (empty perfgraph, no error).  Always clean here.
cd miniapps/seas
make clean
make seas_spatial_dyn_driver
```

## Runtime usage

The p2/O3 sbatch (`jobs/tpv31_spatial/tpv31_p2_aderO3_noflux_50m_flex.sbatch`)
already wires this behind `TPV31_PERFGRAPH=1`: set that env var and it exports a
`CALI_CONFIG` combining a text region tree + a Hatchet JSON. Manually:

```bash
# Text region tree to stderr (.err) at exit:
export CALI_CONFIG="runtime-report(calc.inclusive=true)"

# Call-graph JSON for the Hatchet Python package:
export CALI_CONFIG="hatchet-region-profile(output=${OUT}/caliper_perfgraph.json)"

# Both at once (what the sbatch sets):
export CALI_CONFIG="runtime-report(calc.inclusive=true),hatchet-region-profile(output=${OUT}/caliper_perfgraph.json)"
```

Render the call-graph from the JSON:

```python
import hatchet as ht
gf = ht.GraphFrame.from_caliperreader("caliper_perfgraph.json")
print(gf.tree())        # ASCII call tree with inclusive/exclusive time
# gf.to_dot(...) for Graphviz, or use Thicket for multi-run comparison
```

## Sanity check (catches the stale-object trap)

After a short run with `runtime-report`, the ADER wave + friction regions must
appear. Use regions that actually fire on the ADER path (TPV31) — NOT `Mult`,
which is RK-only and is correctly absent here:

```bash
grep -q 'seas::WaveOperator::ComputeADERFaceFluxRHS' runtime_report.txt || echo "MISSING ADER flux regions — did you 'make clean' (step 3)?"
grep -q 'seas::spatial_dyn::friction_substep'        runtime_report.txt || echo "MISSING friction region"
```

## Notes

- `runtime-report` adds small per-region overhead; for timing-sensitive
  comparison use `hatchet-region-profile` (lower overhead) and compare against a
  no-`CALI_CONFIG` wall-clock baseline.
- The `TPV31_PERFGRAPH` hook currently lives only in the p2/O3 sbatch on the
  `system/spatial_dyn_driver` working tree. Do NOT duplicate it; merge this
  branch into `system/spatial_dyn_driver` so the annotations and the hook reunite
  (see `caliper_perfgraph_plan.md` §Phase 3 "Merge target").
