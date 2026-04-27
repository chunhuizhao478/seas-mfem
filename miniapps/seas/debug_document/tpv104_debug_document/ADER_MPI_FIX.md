# ADER MPI Review (R-1500..R-1510) — Fix Report

**Date:** 2026-04-26
**Companion review:** `ADER_MPI_REVIEW.md` (verdict PASS WITH FIXES; 0 CRITICAL + 4 MODERATE + 6 LOW).
**Scope:** Address actionable findings from the review.  R-1500 is withdrawn; R-1507/R-1509 are declined per review.  R-1504 is deferred (review explicitly says "non-trivial refactor; provide microbenchmark first").

---

## Files modified

| File | Lines (delta) | Findings addressed |
|---|---:|---|
| `dynamic/wave_operator.hpp`            | +29/−2  | R-1501 (mutable buffer members), R-1510 (static-mesh @warning) |
| `dynamic/wave_operator.inl`            | +130/−25 | R-1501 (5 stack→member), R-1503 (cross-ref comments), R-1505 (collective `has_bulk_bg_` check ×2), R-1506 (Q_new sizing hoist), R-1508 (ASSERT→VERIFY ×2), R-1510 (ctor exchange comment) |
| `tests/unit/test_ader_interior_vs_shared_branch_live.cpp` | +110/−10 | R-1502 (parametric ADER-O ∈ {2,3,4} sweep + P=3 fixture + in-element gradient) |

Total: **+269/−37** across 3 files.

---

## Finding-by-finding resolution

### R-1500 — shared-fault one-sided flux (WITHDRAWN by reviewer)
No code change.

### R-1501 MODERATE — per-step heap allocation (LANDED)

**Header.** Added 8 `mutable Vector` members:
- `ader_I_buf_`, `ader_rhs_buf_` (consumed by `AdvanceADER`).
- `ck_D_curr_buf_`, `ck_D_next_buf_`, `ck_dQ_dxd_buf_` (consumed by `ComputeADERTimeIntegrated`).
- `ck_substep_D_curr_buf_`, `ck_substep_D_next_buf_`, `ck_substep_dQ_dxd_buf_` (consumed by `ComputeADERSubStepStates`).

The two CK predictors (TimeIntegrated vs SubStep) get separate buffers so a future driver path that calls both in sequence within one event doesn't fight for the ping-pong slots.

**Implementation pattern.** First call sizes once via `if (size != N) SetSize(N);`; subsequent calls hit the no-op fast path. `D_curr = Q` overwrites prior contents before any read; `D_next = 0.0` zeros the other ping-pong slot before its first accumulation; `dQ_dxd` is fully overwritten by `ApplySpatialDerivative`. No stale-read hazards.

**Correctness.** `mfem::Swap(D_curr, D_next)` swaps internal data pointers between two member Vectors. After swap the members still own valid buffers of size N; on the next call the SetSize check passes, `D_curr = Q` overwrites everything, so the swap state is harmless across calls.

### R-1502 MODERATE — ADER-O sweep on parallel test (LANDED + STRENGTHENED)

**What the review asked for.** Convert `test_ader_interior_vs_shared_branch_live.cpp:main` into a templated/parametric helper and sweep over O ∈ {2, 3, 4}.

**Implementation.** Pulled the fixture-build + serial-vs-parallel slope-extract into `RunAderEquivalenceSlice(rank, nprocs, ader_order)` and call it three times in `main` from the same MPI program.

**Strengthening (beyond the review's literal ask).** The original fixture used per-side **constant** Q via `AddSideConstant`.  Under FE polynomial order P=1 with constant in-element Q, every spatial derivative collapses to zero and `D(k≥1) ≡ 0` numerically — meaning a wrong `(k+2)` factorial coefficient in `ComputeADERTimeIntegrated`'s recursion would still produce zero × wrong = zero on the test's slope output (silent pass).

To make the test coefficient-sensitive at higher orders, two changes were made together:

1. **In-element gradient.** Added `AddInElementGradient` that adds `amp · (j+1)` per local DOF index `j` of each element, so Q is non-constant within each tet.  Applied to VX, SXY, SXZ at 1% of `kVelAmp` / `kStressAmp`.
2. **FE polynomial order P=1 → P=3.**  At P=1 the L2 basis can carry at most a linear field, so `ApplySpatialDerivative` on a piecewise-linear Q produces a piecewise-constant output, and a second derivative collapses to zero — `D(2) ≡ 0` regardless of gradient amplitude.  P=3 lets the recursion produce non-zero `D(2)` and `D(3)`, exercising the `dt³/6` and `dt⁴/24` coefficients at O=3, O=4.  A 2-tet P=3 fixture has 40 DOFs total — comfortably tractable.

**Verification at all three orders.** `mpirun -np 2 ./seas_test_ader_interior_vs_shared_branch_live`:

| ADER-O | Distinct slope from O-1? | Worst rel mismatch (serial vs parallel) |
|---|---|---:|
| 2 | (baseline) | 0.000e+00 |
| 3 | YES (e.g. comp 1: 8.020314e+06 vs 8.020313e+06) | 0.000e+00 |
| 4 | tiny (lower-order CK truncation already small) | 0.000e+00 |

The slope distinction between O=2 and O=3 confirms the higher-order CK contribution is now load-bearing in the test.

### R-1503 MODERATE — CK factorial denominator cross-ref (LANDED)

Added an explanatory comment block at each site (`wave_operator.inl:1019` and `:1108` in the original line numbering, now around `:1018-1037` and `:1124-1143` after the comment additions) explaining:

- `ComputeADERTimeIntegrated` uses `fac *= dt / (k+2)` because `fac` was initialised to `dt` BEFORE the loop (k=0 contribution pre-added with `dt` folded in).
- `ComputeADERSubStepStates` uses `fac[o] *= τ / (k+1)` because `fac[o]` was initialised to `1.0` (k=0 contribution `(τ⁰/0!) · D(0) = D(0)` pre-added unscaled).

The comments cross-reference each other so a future "fix" attempt that "harmonises" the two denominators will see the warning before doing harm.

### R-1504 MODERATE — vector ghost exchange refactor (DEFERRED per review)

Review explicitly says: *"Caveat: shared with `ComputeSharedFaceFluxRHS` which has a deep-copy correctness guard at line 2636-2639 — the refactor must preserve that invariant.  Treat as a separate ticket; provide microbenchmark first."*  Honored.  No code change in this commit.

### R-1505 LOW — collective `MFEM_VERIFY(has_bulk_bg_)` (LANDED)

Wrapped the existing rank-local `MFEM_VERIFY` at both ADER corrector entry points (`ComputeADERFaceFluxRHS` interior, `ComputeADERSharedFaceFluxRHS` shared) with the `SetMixedFluxMode` collective-consensus pattern.  On a parallel mesh:

```cpp
auto &pmesh_consensus = static_cast<const ParMesh &>(mesh_);
const int my_has = has_bulk_bg_ ? 1 : 0;
int min_has = 0;
MPI_Allreduce(&my_has, &min_has, 1, MPI_INT, MPI_MIN, pmesh_consensus.GetComm());
MFEM_VERIFY(min_has == 1, "...");
```

Catches the "rank A skipped `SetAbsorbingBackground` while ranks B+ called it" path before either rank reaches the shared corrector's `ExchangeFaceNbrData`, so the abort is collective rather than asymmetric (which would deadlock the survivors).

Cost: two extra `MPI_Allreduce(MIN)` per `AdvanceADER` call (one per corrector, ~tens of µs each).  Negligible vs. the 12 000 macro-step production scale.  The serial branch keeps the original `MFEM_VERIFY(has_bulk_bg_)` form.

### R-1506 LOW — `Q_new.SetSize` hoist (LANDED)

Moved `Q_new.SetSize(NUM_STATE * ndof_total_)` from the bottom of `AdvanceADER` (post-corrector) to immediately after the aliasing guard.  The bottom site is now annotated as already-sized.  No behavioural change at the current call sites (PML branch writes `rhs`, not `Q_new`) but a future variant that pre-stages into `Q_new` can rely on the contract.

### R-1507 LOW — predictor caching (DECLINED per review)

28 GB cost at 3M tets P=3.  No code change.

### R-1508 LOW — `MFEM_ASSERT` → `MFEM_VERIFY` (LANDED)

Promoted the mixed-element guards at both ADER ghost-layer sites (`ApplySpatialDerivative` ~L860 and `ComputeADERSharedFaceFluxRHS` ~L4235) from `MFEM_ASSERT` (Debug-only) to `MFEM_VERIFY` (always on).  Updated the error messages to be self-documenting.

### R-1509 LOW — redundant size validation (DECLINED per review)

Harmless redundancy.  No code change.

### R-1510 LOW — static-mesh assumption documentation (LANDED)

Added two doc blocks:

1. `wave_operator.hpp` `AdvanceADER`'s docstring: a `@warning R-1510` paragraph explaining that both `ParMesh::ExchangeFaceNbrData()` and `ParFiniteElementSpace::ExchangeFaceNbrData()` are called once at construction, so an AMR / moving-mesh extension MUST re-call them before any subsequent `AdvanceADER` if the partition or face-neighbour graph mutates.
2. `wave_operator.inl:137-148` (the constructor exchange site): a comment explaining the static-mesh contract, the required call order (`ParMesh` first, then `ParFiniteElementSpace`), and pointing to the `AdvanceADER` warning.

---

## Test results

All run from `miniapps/seas/`:

| Test binary | Tests | Result |
|---|---:|---|
| `seas_test_tpv104_substep_iterator`        | 40 | ✅ 40/40 |
| `seas_test_tpv104_substep_iterator_parity` | 22 | ✅ 22/22 |
| `seas_test_tpv104_substep_predictor`       |  6 | ✅ 6/6   |
| `seas_test_tpv104_substep_one_shot_parity` |  4 | ✅ 4/4   |
| `seas_test_tpv104_substep_dispatch_parity` |  1 | ✅ 1/1   |
| `seas_test_tpv104_normal_sign`             | 17 | ✅ 17/17 |
| `seas_test_tpv104_setup`                   | 38 | ✅ 38/38 |
| `seas_test_tpv104_smoke`                   | 53 | ✅ 53/53 |
| `seas_test_tpv104_sigma_n_invariance`      | 20 | ✅ 20/20 |
| `seas_test_tpv104_freeze_sigma_n_gate`     |  6 | ✅ 6/6   |
| `seas_test_tpv104_probe_format`            | 13 | ✅ 13/13 |
| `seas_test_fault_face_flux_ader_equivalence`| 11 | ✅ 11/11 |
| **np=2** `seas_test_ader_interior_vs_shared_branch_live` (O sweep) | 3 (one per ADER order) | ✅ 3/3 (rel=0 at every order, every component) |

`seas_tpv104_driver` builds clean.

---

## Items not addressed

- **R-1504** — vector ghost exchange refactor.  Review marked as deferred ticket.  Expected ~32 s wall-clock saving over a 60 s production run; recommend a microbenchmark before refactoring.  The R-1504 caveat about the deep-copy correctness guard at `ComputeSharedFaceFluxRHS:2636-2639` must be respected by any subsequent refactor.
- **R-1500** — withdrawn by reviewer.
- **R-1507** — declined optimisation (28 GB cost at production scale).
- **R-1509** — declined (harmless redundancy).

---

## Cross-references for follow-on review

* CK factorial denominator comments:
  - `wave_operator.inl:1018-1037` (TimeIntegrated, `(k+2)` denom)
  - `wave_operator.inl:1124-1143` (SubStep, `(k+1)` denom)
* Mutable scratch buffers:
  - declarations: `wave_operator.hpp:599-616`
  - AdvanceADER use: `wave_operator.inl:4523-4533`
  - TimeIntegrated use: `wave_operator.inl:1009-1024`
  - SubStepStates use: `wave_operator.inl:1124-1141`
* Collective `has_bulk_bg_` check: `wave_operator.inl` at the two ADER corrector sites (interior near line 3155, shared near line 4125).
* Q_new sizing hoist: `wave_operator.inl` near AdvanceADER's top, after the aliasing `MFEM_VERIFY`.
* Static-mesh docstring: `wave_operator.hpp:399-413` (AdvanceADER) and `wave_operator.inl:137-148` (ctor).
* Test coverage: `tests/unit/test_ader_interior_vs_shared_branch_live.cpp` — `RunAderEquivalenceSlice` helper + sweep over O ∈ {2, 3, 4} in `main`.
