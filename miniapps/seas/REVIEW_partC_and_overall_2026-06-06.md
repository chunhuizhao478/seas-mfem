# Code Review: Part C (TPV6/7 setup) + overall A+B+C pass

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_unified_bimaterial_volume_and_fault_2026-06-06.md`
- C2 code reviewed: `drivers/spatial_dyn_driver.cpp` (verify-dispatch per-side diagnostic;
  B1 AssignFault call; B3 halfspace build branch), `tests/unit/test_tpv_config_parse.cpp`
  (T_TPV6), `tpv6/configs/tpv6.toml` + `tpv7/configs/tpv7.toml`, `tpv6/mesh/*.geo`,
  the 4 sbatch, `tpv6/visualize_results.py`, `tpv6/benchmark_data/README.md`.
- Reference: `tpv6/benchmark_document/2007RuthRalphletter2.pdf` (spec) +
  `tpv6/benchmark_data/scec_drdg3d` (DRDG3D reference traces).
- Verified CORRECT (not findings):
  - **near/far material (spec p.3 + drdg3d data):** near=STRONG (vp2=6000), far=WEAK
    (vp1=3750 TPV6 / 5000 TPV7).  `tpv6.toml`/`tpv7.toml` set vp_near=6000, vp_far per
    problem — MATERIAL ASSIGNMENT CORRECT.  Empirically the far (weak) side moves ~3-4x
    faster in the reference; `--verify-dispatch` prints Zp_minus(near)=1.602e7 != Zp_plus(far)=8.34e6.
  - Configs spec-match: tau_strike 70/81.6 MPa, sigma_n 120, mu_s 0.677/mu_d 0.525/d_c 0.40,
    3000 m square patch, barriers (|x|>15km, z<-15km), free surface; TPV7 far triple.
  - Nucleation = TPV205 static stress patch (no [nucleation] block); the abandoned
    `instantaneous_overstress_square` kind was fully REVERTED (compiles clean, no leftovers).
  - B3 driver branch lifetime: `halfspace_wrapper` is declared with `depth_profile_wrapper`
    (outlives `wave_ptr`/`material`); the MaterialField borrows its FunctionCoefficients.
  - Byte-exact contract intact: homog_equivalence 18/18; B1 no-op on scalar + Zp-equal on
    depth-symmetric; AssignFault + the diagnostic gated so TPV31/TPV102/TPV205/BP5 unchanged.
  - Both run arms construct on the 200 m mesh (`--dry-run`); config-parse 140/140.

## Findings

### [R-001] MODERATE — spatial_dyn_driver.cpp: the verify-dispatch per-side fault banner counts rank-0-LOCAL fault DOFs, not global

**Category:** BUG (diagnostic correctness under MPI)

**Description:**
The new `[verify-dispatch] fault per-side` line scans `dof_data` inside
`if (verify_dispatch && rank == 0)`.  `dof_data` holds only THIS rank's local fault
DOFs, so the printed `N/M DOFs` count and the `Zp_plus/Zp_minus` min/max reflect
rank 0's partition only.  At np>1 the fault may be unevenly partitioned (or absent
from rank 0 -> banner silently skipped), so the headline acceptance check
("banner shows Zp_plus != Zp_minus") is unreliable for a parallel dry-run.

**Trigger:** `mpirun -np >1 ... --dry-run --verify-dispatch` where rank 0 owns few/no fault DOFs.

**Actual behavior:** rank-0-local count + min/max (correct only at np=1).

**Expected behavior:** global count + global Zp min/max (or the line explicitly labelled "rank-0 local").

**Suggested fix:** compute the local tallies on ALL ranks (verify_dispatch is
rank-uniform), MPI_Reduce to rank 0, then print on rank 0:
```cpp
// (outside the rank==0 block; all ranks participate)
long long n_bimat_l = 0, n_dof_l = (long long)dof_data.size();
real_t zpmn_l =  std::numeric_limits<real_t>::max(), zpmx_l = 0.0;
real_t zmmn_l =  std::numeric_limits<real_t>::max(), zmmx_l = 0.0;
for (const auto &d : dof_data) { /* accumulate as today */ }
long long n_bimat = n_bimat_l, n_dof = n_dof_l;
real_t zpmn = zpmn_l, zpmx = zpmx_l, zmmn = zmmn_l, zmmx = zmmx_l;
#ifdef MFEM_USE_MPI
MPI_Reduce(&n_bimat_l, &n_bimat, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
MPI_Reduce(&n_dof_l,   &n_dof,   1, MPI_LONG_LONG, MPI_SUM, 0, comm);
MPI_Reduce(&zpmn_l,&zpmn,1,MPITYPE,MPI_MIN,0,comm); MPI_Reduce(&zpmx_l,&zpmx,1,MPITYPE,MPI_MAX,0,comm);
MPI_Reduce(&zmmn_l,&zmmn,1,MPITYPE,MPI_MIN,0,comm); MPI_Reduce(&zmmx_l,&zmmx,1,MPITYPE,MPI_MAX,0,comm);
#endif
if (verify_dispatch && rank == 0 && n_dof > 0) { /* print n_bimat/n_dof + ranges */ }
```
(MPITYPE = MFEM's `MPI_DOUBLE`/real_t type.)

**Test case:** not unit-testable without an np>1 fixture; verified by reasoning + the
np=1 dry-run (78387/78387).  MODERATE because it is the acceptance-criterion banner.

---

### [R-002] LOW — tpv6/visualize_results.py: only the strike (h-*) field family is auto-overlaid; reference n-stress sign not flipped

**Category:** QUALITY

**Description:**
The viz plots one `--field` at a time (good) and reads the drdg3d reference, but when
overlaying THIS run's output it does not flip the drdg3d compression-NEGATIVE n-stress
to this repo's compression-POSITIVE convention, so an `--field n-stress` overlay would
show the reference at -120 vs the run at +120.

**Suggested fix:** when `--field n-stress`, negate the reference column before plotting
(documented in README); or annotate the axis. LOW — the writer that produces the run
side is itself deferred.

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — rank-0-local verify-dispatch banner)
- Low issues: 1 (R-002 — viz n-stress sign on overlay)
- Plan compliance: FULL for the verifiable C2 core (configs/meshes/jobs/dispatch);
  one DOCUMENTED plan correction (contrast guard reclassifies a few welded-interface
  faces, not 0) and one DOCUMENTED deferral (the per-side station writer).
- Verdict: PASS WITH FIXES — no critical bugs.  Apply R-001 (acceptance banner
  correctness under MPI); R-002 is cosmetic.

## Overall (A+B+C) — holistic check
- **Byte-exact contract HOLDS** for every shipped problem: Part A contrast guard
  defaults OFF (tol<0); Part B per-side fault is no-op on scalar + Zp-equal on
  depth-symmetric faults (TPV31/TPV102/TPV205/BP5); Part C adds only new TPV6/7
  artifacts + a flag-gated diagnostic.  Evidence: homog_equivalence 18/18,
  contrast_guard 16/16, config-parse 140/140 (existing configs unchanged).
- **Composition validated end-to-end:** `--dry-run` on the real 200 m bi-material
  mesh constructs both arms (matrix + per-side fault + halfspace material), confirming
  A (guard), B (per-side Riemann + B1 assignment + B2 relaxation + B3 material) and C
  (config/mesh/jobs) compose.
- **No reverted-fix violations** (CLAUDE.md): the only revert was the just-added,
  unused square nucleation kind (per the user's "don't duplicate code"); no prior
  fix was touched.
- **Remaining work (documented):** the per-side ON-FAULT station writer
  `dynamic/tpv6_stations.hpp` (per-side disp+vel + traction; strong side="nearside";
  flip n-stress sign for the drdg3d overlay) + its driver tag case at
  spatial_dyn_driver.cpp:3110; and the Frontera handedness validation (flip
  `[material.halfspace_across_fault].n_y` if the rupture is strike-mirrored vs drdg3d).

## Unreviewed Areas
- Full TPV6/7 dynamic rupture vs drdg3d — Frontera (no local full-mesh sim).
- np>1 partition behavior of the fault per-side assignment (shared-fault R-101 deferred).
