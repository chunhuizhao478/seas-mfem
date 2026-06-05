# Code Review + Fix Report: Phase-0 fixes & Phase-1 fault-flux over-integration (2026-06-02)

## Scope
- Plan: `seisol_overintegration_resample_speckle_2026-06-02.md` §6 **Phase 1** (over-integration of the fault flux — the prerequisite for resample).
- Prior review consumed: `REVIEW_phase0_code.md` (findings R-001…R-005).
- Files changed:
  - `dynamic/wave_operator.hpp` (+39) — `SetFaultOverint`/`GetFaultOverint`/`FaultFaceQuadDegree`, member `fault_overint_k_`, helper `RebuildFaultQuadrature_`.
  - `dynamic/wave_operator.inl` (+189/−56) — ctor factored to `RebuildFaultQuadrature_()`; new helper + setter; 8 fault-quadrature sites made fault-aware.
  - `tests/unit/test_fault_planar_serial.cpp` — R-001…R-005 fixes + `--fault-overint` wired live.
  - `tests/unit/test_fault_overint.cpp` (new) — Phase-1 unit test (tiny 2-tet fixture).
  - `Makefile` — `seas_test_fault_overint` target + `test-fault-planar-serial-fast` (R-005).
- Build/test environment per `[[reference_worktree_build_mfem_override]]`: `MFEM_DIR=../.. MFEM_BUILD_DIR/INC_DIR/LIB_DIR=<main checkout>`, conda `mfem-dev`.

---

## Part A — Phase-0 review fixes (R-001…R-005)

| ID | Severity | Resolution | Verified |
|----|----------|-----------|----------|
| R-001 | MODERATE→**corrected** | dt auto-scales with `he` and the `3/(2N+1)` order de-rating when `--dt` not given; CFL backstop `MFEM_VERIFY(o.dt ≤ 2·ComputeMaxDt(1))`. **Deviation from the reviewer's suggested `ComputeMaxDt(0.5)` guard** — see note below. | `--dt 0.2`→aborts before stepping (exit 1); default & `--he 250` run; `ComputeMaxDt(1)=0.069 s` |
| R-002 | LOW | `--nsteps` help + header example aligned to the code default (8000). | `--help` shows 8000 |
| R-003 | LOW | `--ny-half N` knob (default **2**); `BuildFaultSlab` generalized to `2·ny_half` layers, fault at `y=ny_half·he`, symmetric = lower-half `+1`/upper-half `−1`; `VerifyMirrorSymmetry` uses `Ly=2·ny_half·he`. | mirror symmetry empirically confirmed for `ny_half∈{1,2,3}` (independent reviewer) |
| R-004 | LOW | `#include <array>` added. | compiles |
| R-005 | LOW | `test-fault-planar-serial-fast` CI smoke; documented that the full run is not in the `test` aggregate. | n/a |

**R-001 deviation (the reviewer's own suggested fix was wrong).** The Phase-0 review proposed `MFEM_VERIFY(o.dt ≤ wave.ComputeMaxDt(0.5))`. Measured at runtime, `h_min` is the tet **inscribed diameter** `6·Vol/A_total ≈ 0.41·he ≈ 414 m`, so `ComputeMaxDt(0.5)=0.5·h_min/cp ≈ 3.45e-5 s` — **below** the validated default `dt=1e-4`. That guard would have aborted the byte-exact default run (which sits at cfl ≈ 1.5e-3 of `ComputeMaxDt(1)`, i.e. deeply sub-CFL). The original review's "`--he 250` blows up" trigger was therefore a unit-arithmetic error: `dt=1e-4` has ~600× CFL margin, so `--he` would have to drop to ~1 m to destabilize. The replacement (auto-scale + a generous `2×ComputeMaxDt(1)` backstop) keeps the default byte-exact, makes non-default `--he`/`--order` runs stay CFL-stable, and still catches a grossly over-CFL explicit `--dt` (verified: `--dt 0.2` aborts; `--dt 1e-4` default does not). The order de-rating also matches the plan's planned p1-vs-p2 acceptance.

---

## Part B — Phase-1 implementation (fault-flux over-integration)

**Design.** A single knob `WaveOperator::SetFaultOverint(k)` (default off, `k=0`) sets `fault_overint_k_` and rebuilds the fault quadrature. `FaultFaceQuadDegree() = 2*(order_+k)` is the fault-face exactness degree; at `k=0` it is exactly `2*order_`, so **every** fault-quadrature site is byte-identical to pre-Phase-1 *by construction*. The 8 fault-quadrature sites were split per an exhaustive classification:
- **FAULT-ONLY (unconditional `FaultFaceQuadDegree()`):** ctor `RebuildFaultQuadrature_` (`nbf_per_face_` + `FaultBasis` QP data), `EvaluateBulkAtFaultQPsCanonical` (interior+shared), `ExchangeAndPairSharedFaultQPs`, `VerifySharedFaultDOFDataConsistency`.
- **MIXED (fault-aware ternary `face_is_fault/sf_fault ? FaultFaceQuadDegree() : 2*order_`):** `ComputeFaceFluxRHS`, `ComputeADERFaceFluxRHS`, `ComputeSharedFaceFluxRHS`, `ComputeADERSharedFaceFluxRHS`.
- **Untouched (element-volume `fe->GetGeomType()`):** `ComputeVolumeRHS`, `ApplySpatialDerivative`, `AssembleElementMassInverse`, `ApplyPMLDamping`, the deriv/volume-op builders.

Friction parameters are assigned per over-integration GP automatically: the harness/driver build `fault_coords` + `DOFData` from `FaultFaceQuadDegree()`, so `InitializeFaultDOFs` (which reads per-QP coords) lands one parameter set per GP.

**Guards added.** `SetFaultOverint` aborts if called after `SetFaultDOFData` (would mis-size DOFData), and rejects `k>0` combined with mixed-flux / precomputed-face-flux (those cache at `2*order_`; out of Phase-1 scope).

### Adversarial review — CRITICAL bugs: **none found**

Independently re-derived (two reviewers: author + an independent agent over the full diff). Verdicts by risk category:

- **(a) Byte-exactness at k=0 — CLEAR.** All 14 `IntRules.Get` sites accounted for; no element-volume site changed; no fault site missed. At `k=0` the ternary collapses (`FaultFaceQuadDegree()==2*order_`) so both branches select the same rule.
- **(b) Predicate divergence — CLEAR.** At each MIXED routine the rule-degree predicate is *identical* to the per-QP handling predicate (`ComputeFaceFluxRHS`/`ComputeADERFaceFluxRHS`: `face_is_fault ≡ is_fault` because the routine's `is_fault` is only reached under `!is_boundary`; the two shared routines reuse the *same* `sf_fault` variable). No path gives over-integrated QP count with bulk handling, so `fault_face_dof_offset_[f]+q` / `shared_fault_dof_offset_[sf]+q` can't go out of range.
- **(c) `RebuildFaultQuadrature_` == original ctor block — CLEAR.** Verbatim move with the single `2*order_ → FaultFaceQuadDegree()` change; same `ref_normal=(0,-1,0)`, `up=(0,0,1)`, same `face_geom` logic, same `ComputeQPBasis`/`ComputeQPBasisShared(...)`. `FaultBasis::ComputeQPBasis*` use `qp_data.resize(nq)` (overwrite, not append) ⇒ ctor build and any `SetFaultOverint(0)` rebuild are bit-identical.
- **(d) Scope/lifetime — CLEAR.** No use of the removed `face_ir`/`face_geom` locals survives; `fault_basis_` deref is guarded; helper runs after `make_unique<FaultBasis>()`.
- **(e) dt auto-scaling — CLEAR.** Factor is exactly `1.0` at default order=1/he=1000 (byte-exact); applied once in `main()` before Part B copies `o`; backstop cannot false-abort the default.
- **(f) `ny_half` generalization — CLEAR.** Mirror symmetry holds for the symmetric mesh and fails for the asym mesh at `ny_half∈{1,2,3}` (empirically reproduced); `nyv=2·ny_half+1`, fault at `0.5·Ly`, `VerifyMirrorSymmetry` uses matching `Ly`.
- **(g) Unit-test OOB/uninit — CLEAR.** All buffers sized from the mutually-consistent `GetNbfPerFace()`/`FaultFaceQuadDegree()`; `SetFaultOverint` precedes `SetFaultDOFData`.

**One non-bug observation:** on a rank with no fault faces, `SetFaultOverint(k)` leaves `nbf_per_face_==0` while `FaultFaceQuadDegree()` reports `2*(order+k)`. Harmless (no QPs to size; matches pre-existing ctor behavior).

### Validation evidence
`make test-fault-overint` → **10/10 pass** on the 2-tet fixture, including:
- **T5 (byte-exact-when-off, runtime):** one `AdvanceADER` step with `SetFaultOverint(0)` is **bit-identical** (`max|dQ| == 0`) to never calling it.
- **T6 (knob live):** `k=1` changes the radiated flux; `GetNbfPerFace()` grows 3→6 (TRIANGLE degree 2→4); `FaultBasis` qp_data resized accordingly.

---

## Verdict
**PASS.** No critical bugs. The k=0 / "knob-off" regression contract for TPV102/104/205/BP5/spatial is preserved **by construction** (every fault-quad site collapses to `IntRules.Get(geom, 2*order_)` at k=0) and is runtime-confirmed on a fault fixture (T5). The live over-integration plumbing is internally consistent across all 8 fault sites + the harness.

## Residual validation (deferred — cannot run locally)
Per `[[feedback_no_local_reproducer]]`, the production reproducer was NOT run locally. The following need Frontera (plan §7) or an approved run:
1. **Production byte-exact-off** on TPV102/104/205/BP5/spatial — guaranteed by construction + T5, but not run end-to-end on the production driver paths.
2. **MPI (np>1) over-integration** — the shared-fault inline path (`ComputeADERSharedFaceFluxRHS`, `ExchangeAndPairSharedFaultQPs`) is over-int-aware and consistent by construction (same `k`→same rule→same QP ordering on all ranks), but not run at np>1.
3. **The actual physics payoff** — that `--fault-overint` reduces the on-fault σ_n speckle/drift on the Phase-0 harness, and the predicted TPV102/TPV31 flattening. This is the Phase-0 A/B the harness exists to run.
4. **Substep-iterator path (Phase 3 dependency):** `tpv*_substep_iterator.cpp` still build their own fault QP arrays at `2*order` and are NOT over-int-aware. This is safe today (over-int is off in production; the harness uses the inline `AdvanceADER` path), but Phase 3 must thread `FaultFaceQuadDegree()` through the iterators before over-integration can drive the production substep dispatch.
5. **`ny_half=2` default** changes the harness's default mesh; the Part-B mirror assertions self-check at runtime, but the full acceptance run (does the drift still reproduce at 2-thick?) is unverified.

## Low / future items (non-blocking)
- `SetMixedFluxMode` / `UsePrecomputedFaceFluxes` have no reciprocal guard against an already-set `fault_overint_k_>0` (the guard lives only in `SetFaultOverint`). Add a symmetric check when those paths are made over-int-aware.
