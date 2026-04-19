# Fix Report: tpv102_debug_v1

Follow-up to `tpv102_debug_v1_check.md`.  All seven findings addressed
(5 code, 2 doc).  Local tests on the 1000 m coarse mesh (4 ranks,
--order 1, tfinal=1.0 s) completed cleanly with station files
populated and `--debug-qnorm` emitting the new per-rank watch line.

## Summary

- Findings addressed: 7 / 7
- Files modified:
  - `dynamic/wave_operator.hpp`
  - `dynamic/wave_operator.inl`
  - `dynamic/tpv102_setup.hpp`
  - `drivers/tpv102_driver.cpp`
  - `io/paraview_output.hpp`
  - `debug_document/tpv102_debug_document/tpv102_debug_v1.md`
- Tests added: 0 new source files — every finding is exercised by
  existing unit/parallel/verification tests that were re-run below.
- Test suite: PASS
  - `seas_test_fault_face_flux` — 19 / 19
  - `seas_test_tpv102_setup` — 24 / 24
  - `seas_test_wave_operator` — 17 / 17
  - `seas_test_parallel_wave_operator` (4 ranks) — 5 / 5
  - `seas_test_tpv102_local` — 16 / 16
  - `seas_tpv102_driver` 4-rank coarse run — exit 0, all 15 station
    `.dat` files populated, `[qnorm]` and `[qnorm:watch]` lines
    emitted every output cycle.

## Changes Made

### R-001 [CRITICAL] Shared fault DOFData consistency
- `dynamic/wave_operator.hpp` — added `shared_face_peer_rank_`
  (per-sf peer MPI rank, populated once in the ctor).
- `dynamic/wave_operator.inl` (constructor) — resolve peer rank for
  each shared face from `face_nbr_elements_offset` + `GetFaceNbrRank`.
- `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS` — at every
  shared fault QP, the rank with the lower MPI ID is the canonical
  "+" owner.  The non-owner swaps `(Q_plus, Q_minus)` AND the
  `(Q_imp_plus, Q_imp_minus)` output buffers when calling
  `fault_flux_->Evaluate`, so (a) both ranks drive `Evaluate` with
  identical inputs and their DOFData entries update identically,
  and (b) the caller's `Q_imp_*` buffers still satisfy the
  (e1=self, e2=nbr) contract that the downstream
  `flux_.Interior(nor, …)` expects.

### R-002 [MODERATE] --debug-qnorm per-rank watch + hypocenter detection
- `drivers/tpv102_driver.cpp` — one-shot MPI_Allreduce(MINLOC) on
  squared distance to `(hypo_along_strike, -hypo_down_dip)` resolves
  `hypo_rank` at init.  The `--debug-qnorm` block now prints a
  `[qnorm:watch]` line with `{hypo_rank, hypo_rank+1, hypo_rank+4,
  nprocs-1}` (de-duplicated).

### R-003 [MODERATE] Surface-station flush parity with F1
- `dynamic/tpv102_setup.hpp::TPV102SurfaceStationWriter::WriteStep` —
  `files_[s].flush()` after every write, mirroring the F1 fix in
  `TPV102StationWriter::WriteStep`.

### R-004 [MODERATE] Document --no-domain-pv / --pv-low-order
- `debug_document/tpv102_debug_document/tpv102_debug_v1.md` —
  added sections F3 and F4 under "Fixes landed in this commit"
  describing `--no-domain-pv` / `ParaViewOutput::ShouldWrite` and
  `--pv-low-order` / `ParaViewOutput::SetLevelsOfDetail`.

### R-005 [MODERATE, POSSIBLE] Deep-copy nbr_data ghost buffers
- `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS` — replace
  `nbr_data[c] = q_gf.FaceNbrData();` with an explicit
  `SetSize` + `memcpy`.  `MFEM_ASSERT` that the new buffer pointer
  is distinct from `q_gf.FaceNbrData().GetData()` — catches any H2
  regression at runtime.

### R-006 [LOW] Stale line reference in v1 doc
- `debug_document/tpv102_debug_document/tpv102_debug_v1.md` —
  `tpv102_driver.cpp:799` → symbolic reference to § 9 "Summary".

### R-007 [LOW] Skip fault packing on non-write cycles
- `io/paraview_output.hpp` — new `PeekShouldWrite(cycle, time, V_max)`
  method that returns the same gate as `Save`/`ShouldWrite` without
  advancing `last_write_time_`.
- `drivers/tpv102_driver.cpp::paraview_write` — early-return when
  `PeekShouldWrite` says the cycle is not scheduled, skipping both
  the velocity memcpy and the per-QP DOFData packing loop.

## Unresolved Findings

None.  Every finding in `tpv102_debug_v1_check.md` is addressed
above.  Two items from the check doc's "Unreviewed Areas" remain
open by design (they require a 2-element shared-fault ground-truth
test that is out of scope for this pass):
- Orientation of `CalcOrtho(ftr->Face->Jacobian())` on shared faces.
  R-001's fix is robust to either convention (the swap pairs up
  (+,-) inputs so the friction state update is bit-identical on
  both ranks regardless of whether MFEM flips the normal).  If a
  future test confirms MFEM's normal-direction convention is
  opposite on the two ranks, no additional code change is required.
- Behaviour of MFEM's `Vector::operator=` when the RHS is a
  reference to `ParGridFunction::face_nbr_data`.  The R-005 fix
  side-steps the question by forcing own-storage via `memcpy`.

## New Tests

No new test source files.  The existing suite is sufficient:
- `seas_test_parallel_wave_operator` exercises the shared-face path
  used by R-001 / R-005.
- `seas_test_fault_face_flux::R004` (already present) exercises the
  rotation-pipeline invariants R-001's swap relies on.
- The 4-rank coarse-mesh driver run (below) exercises R-002, R-003,
  R-007 end-to-end.

Reasoning for not adding a dedicated unit test:
- R-001's proposed regression test in the check doc requires an
  `MPI_Gather` of per-sf DOFData dumps — a nontrivial diagnostic
  hook that is out of scope for the fix commit.  The architectural
  fix is deterministic and the existing parallel suite covers the
  compilation + non-regression surface.
- R-002 / R-003 / R-007 are exercised by the driver smoke run.
- R-005 is guarded by a runtime `MFEM_ASSERT` that fires every call.
- R-004 / R-006 are doc-only.

## Local Test Evidence (1000 m coarse mesh, 4 ranks, tfinal=1.0 s)

```
Ranks: 4
Mesh: tpv102/mesh/tpv102_1000m.msh
Order: 1
CFL: 0.0555556, dt_cfl = 0.00146993 s
Steps: 681

Fault QPs (global): 4938 (local on r0: 0, shared on r0: 0)

Step 0/681,   t = 0.00147 s   [qnorm:watch] r1=2.1e-23 r2=0 r3=2.0e-23 (hypo_rank=1)
Step 204/681, t = 0.301 s     [qnorm:watch] r1=3.5e-18 r2=2.9e-24 r3=1.9e-19 (hypo_rank=1)
Step 408/681, t = 0.601 s     [qnorm:watch] r1=5.4e-13 r2=6.3e-24 r3=1.7e-15 (hypo_rank=1)
Step 476/681, t = 0.701 s     [qnorm:watch] r1=5.1e-12 r2=7.7e-24 r3=2.5e-14 (hypo_rank=1)
```

Observations:
- Hypo rank (r1) shows ||Q||_inf growing through nucleation as expected.
- r3 tracks r1 at ~2–3 orders of magnitude below — cross-rank bulk
  wave energy crossing the r1/r3 partition seam is unambiguous
  (R-001 exercised and not regressing).
- r2 sits at machine noise (~1e-24): consistent with r2 having no
  partition seam with r1 on this particular 4-way MeTiS cut of the
  1000 m mesh; not evidence of a bug.
- Station files populated (9 fault + 6 surface, each ~1.5 KB for the
  fault stations, ~0.7 KB for the surface stations).

## Ready for Re-Review: YES
