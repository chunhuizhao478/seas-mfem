# Known Issue: Antiplane Template Instantiation Broken on This Branch

**Status:** Deferred — to be resolved when the working antiplane branch is merged into `feature/paraview-compaction`.
**Date opened:** 2026-05-11
**Branch:** `feature/paraview-compaction`
**Affected platforms:** All (Frontera, local mfem-dev). First observed on Frontera while building Phase 6 ParaView-compaction.

## Symptom

Any compilation unit that instantiates

```cpp
mfem::seas::SEASQuasiDynamicOperator<MeshType,
                                     mfem::seas::AntiplaneDomainOperator<MeshType>,
                                     mfem::seas::RateStateFaultOperator<MeshType, 1>>::Mult(...)
```

fails with:

```
solver/seas_operator.hpp(342): error:
  class "mfem::seas::AntiplaneDomainOperator<...>"
  has no member "ComputeTractionDiagnostics"

solver/seas_operator.hpp(478): error:
  class "mfem::seas::AntiplaneDomainOperator<...>"
  has no member "IsFirstStepDebugEnabled"
```

## Root cause

`solver/seas_operator.hpp::SEASQuasiDynamicOperator::Mult` calls two
debug methods on `domain_->`:

| seas_operator.hpp line | Call                                                  | Purpose                          |
|------------------------|-------------------------------------------------------|----------------------------------|
| 342                    | `domain_->ComputeTractionDiagnostics(...)`            | BP5 face-tracer split traction   |
| 478                    | `domain_->IsFirstStepDebugEnabled()` (gate)           | First-step debug print           |

Both methods exist on `ElasticityDomainOperator` (BP5 3D path) but were
**never added to `AntiplaneDomainOperator`** (BP1/BP2 2D antiplane path).
The template `Mult` body calls them unconditionally on `domain_->`, so any
instantiation with `AntiplaneDomainOperator` fails at compile time.

Introduced by:
- `8cb7c00` (2026-04-05, "Add accepted-step face tracing for BP5 diagnostics") — added
  `ComputeTractionDiagnostics` and the calling site at line 342.
- `b631c63` (2026-04-06, "Add per-stage traction and V norm diagnostics") — added
  `IsFirstStepDebugEnabled` and the calling site at line 478.

Neither commit touched `domain/antiplane_operator.hpp`, leaving BP1/BP2 build
paths broken from those commits onward.

## Affected build targets (do not build until fixed)

All 14 cpp files that combine `SEASQuasiDynamicOperator` with
`AntiplaneDomainOperator`:

| Source file                                                    | make target                              | Filter pattern   |
|----------------------------------------------------------------|------------------------------------------|------------------|
| `pseas.cpp`                                                    | `seas_pseas`                             | `*pseas*`        |
| `tests/verification/bp1_verification_full.cpp`                 | `seas_bp1_full`                          | `*bp1*`          |
| `tests/verification/bp2_verification_full.cpp`                 | `seas_bp2_full`                          | `*bp2*`          |
| `tests/verification/bp2_verification_first_cycle.cpp`          | `seas_bp2_verify`                        | `*bp2*`          |
| `tests/verification/bp2_benchmark_parallel.cpp`                | `seas_bp2_benchmark_parallel`            | `*bp2*`          |
| `tests/verification/bp2_serial_smoke.cpp`                      | `seas_bp2_serial_smoke`                  | `*bp2*`          |
| `tests/unit/test_bp2_short.cpp`                                | `seas_test_bp2_short`                    | `*bp2*`          |
| `tests/unit/test_quasi_dynamic.cpp`                            | `seas_test_quasi_dynamic`                | explicit         |
| `tests/unit/test_io.cpp`                                       | `seas_test_io`                           | explicit         |
| `tests/unit/test_checkpoint.cpp`                               | `seas_test_checkpoint`                   | explicit         |
| `tests/parallel/test_scaling.cpp`                              | `seas_test_scaling`                      | explicit         |
| `tests/parallel/test_serial_parallel_consistency.cpp`          | `seas_test_serial_parallel_consistency`  | explicit         |
| `tests/parallel/test_parallel_fault.cpp`                       | `seas_test_parallel_fault`               | explicit         |
| `tests/parallel/test_br2_consistency.cpp`                      | `seas_test_br2_consistency`              | explicit         |

`build_frontera.sh` filters these by combining a wildcard pattern
(`*bp1*|*bp2*|*antiplane*|*pseas*`) with an explicit skip list for the
seven targets the pattern doesn't catch.

## Why this is deferred (not fixed in-tree)

User has a working antiplane branch (with `ComputeTractionDiagnostics` /
`IsFirstStepDebugEnabled` stubs or with the template guarded). That branch
will be merged into `feature/paraview-compaction` once the ParaView
compaction work is done. Patching `seas_operator.hpp` here in isolation
risks conflicts with that merge — `seas_operator.hpp` is flagged
"extreme care" in `miniapps/seas/CLAUDE.md` ("Wrong phase order ->
bad equilibrium").

For now, **focus on 3D (BP5 + TPV*)**, which is unaffected:
`ElasticityDomainOperator` has both methods.

## Possible fixes (for when the merge happens)

Listed in order of decreasing invasiveness — the merge from the working
antiplane branch should pick one. Do **not** apply these now without
checking what the antiplane branch already does.

1. **Add no-op stubs to `AntiplaneDomainOperator`** (least invasive,
   least risk). Mirror the `ElasticityDomainOperator` signatures and
   make them no-ops:
   ```cpp
   // domain/antiplane_operator.hpp
   void ComputeTractionDiagnostics(/* same signature as ElasticityDomainOperator */) {
       // Antiplane has no BP5 face-tracer wiring; no-op.
   }
   bool IsFirstStepDebugEnabled() const { return false; }
   ```
   Doesn't change numerics for BP1/BP2; the calls are gated by
   `face_tracer_->IsActive()` and `stage_count < 10 && ... `, both of
   which are false in antiplane runs.

2. **Guard the calls with `if constexpr`** (touches the "extreme care"
   file but only adds compile-time guards):
   ```cpp
   // solver/seas_operator.hpp:339
   if constexpr (detail::has_compute_traction_diagnostics_v<DomainOpType>) {
       if (face_tracer_ && face_tracer_->IsActive()) {
           domain_->ComputeTractionDiagnostics(...);
           ...
       }
   } else {
       // antiplane fallback: skip the diagnostics call entirely
   }
   ```

3. **Specialize `Mult` for `AntiplaneDomainOperator`** (most invasive,
   most code duplication; not recommended).

## How to verify the fix later

When the working antiplane branch is merged, remove the skip from
`build_frontera.sh`:

```bash
# build_frontera.sh — restore make all behavior
make all -j"${JOBS}"
make seas_tpv102_driver seas_tpv104_driver seas_tpv205_driver -j"${JOBS}"
```

Then re-run the Frontera build and confirm all 14 targets compile.
Spot-check `make test-friction-law` and `make test-bp5-smoke` still pass.

## Related references

- `miniapps/seas/CLAUDE.md` § "Files Requiring Extreme Care": flags
  `solver/seas_operator.hpp` as high-risk.
- `domain/elasticity_operator.hpp:221, 310` — the canonical method
  definitions that antiplane must match.
- Commit `8cb7c00` (2026-04-05), `b631c63` (2026-04-06) — introduced the
  asymmetry.
