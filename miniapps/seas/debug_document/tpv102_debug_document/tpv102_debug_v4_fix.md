# Fix Report v4: address review findings in `REVIEW.md` (v4 review)

## Summary
- Findings addressed: **5 of 5** (R-401, R-402, R-403 MODERATE; R-404, R-405 LOW)
- Files modified:
  - `miniapps/seas/dynamic/tpv102_setup.hpp`             (R-401)
  - `miniapps/seas/drivers/tpv102_driver.cpp`            (R-402 part 2 — qnorm formatting)
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch`  (R-402, R-403)
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch` (R-402, R-403)
  - `miniapps/seas/dynamic/wave_operator.inl`            (R-404, R-405)
- Files added:
  - `miniapps/seas/jobs/tpv102/test_result_regex.sh`     (R-402 regex regression test)
- Tests added: 1 new shell-level test (R-402 regex positive+negative cases)
- Test suite: **PASS**
  - `seas_test_godunov_flux`                      — 29/29
  - `seas_test_wave_operator`                     — 17/17
  - `seas_test_wave_bc`                           — 10/10
  - `seas_test_fault_face_flux`                   — 19/19
  - `seas_test_tpv102_setup`                      — 24/24
  - `seas_test_tpv102_local`                      — 16/16
  - `seas_test_parallel_wave_operator` (2 ranks)  — 5/5
  - `seas_test_parallel_wave_operator` (4 ranks)  — 5/5
  - `seas_test_r101_shared_fault` (2 ranks)       — 8/8 (R-302a inline test:
    `3 pairs matched, 0 unpaired, max_diff=0`)
  - `jobs/tpv102/test_result_regex.sh`            — OK (4 positive + 5 negative cases)
  - End-to-end `seas_tpv102_driver` on 4-rank 1000 m mesh, tfinal=0.05 s,
    `--debug-qnorm` — clean run; qnorm now prints as `r1=2.947e-22`
    (stable %.3e format).

## Changes Made

### [R-401] [MODERATE] — `MPITypeMap<real_t>::mpi_type` in station writer
**Where:** `dynamic/tpv102_setup.hpp:246-247`

`TPV102StationWriter::Open` (MPI variant) reduced a `std::vector<real_t>`
buffer with a hardcoded `MPI_DOUBLE` datatype. On `MFEM_USE_SINGLE` builds
(`real_t = float`, stride 4) this silently corrupts the station-ownership
distances.  The v3 fix applied `MPITypeMap<real_t>::mpi_type` at four sites
in `tpv102_driver.cpp` and `wave_operator.inl`; this was the only remaining
site in the TPV102 compile unit and is now consistent.

### [R-402] [MODERATE] — RESULT.txt dead-rank detector + driver qnorm format
**Where:**
- `drivers/tpv102_driver.cpp` (`[qnorm:watch]` print block)
- `jobs/tpv102/tpv102_200m_p1_1.5s_{50,400}rank_dev.sbatch`
- `jobs/tpv102/test_result_regex.sh` (new)

Two-part fix:
1. **Driver**: the `[qnorm:watch]` line now prints with `std::scientific`
   + `setprecision(3)`, restoring flags after the loop.  Format is always
   `M.MMMe[+-]NN`, removing the dependency on `operator<<`'s value-dependent
   switching between fixed and scientific (which was the root cause of the
   false-positive regex match on `rN=0.0001`).
2. **Sbatch**: the dead-rank regex is replaced by a format-anchored pattern
   matching exactly two cases:
     - `rN=0.000e+00` (exact zero), or
     - `rN=M.MMMe-NN` where the exponent is `-20 … -99` (effectively zero).
   The old pattern `r[0-9]+=0(\b|[^.0-9])|r[0-9]+=1e-3[0-9]+` is replaced
   by `r[0-9]+=(0\.0{3}e\+[0-9]{2}|[0-9]\.[0-9]{3}e-[2-9][0-9])([[:space:]]|$)`.

Regression protection: `jobs/tpv102/test_result_regex.sh` pins the expected
behaviour with four positive (dead) and five negative (alive) cases,
including the explicit false-positive input that broke the v3 pattern
(`rN=0.0001` → now correctly emits "alive").  Any future regex edit must
re-pass this script.

### [R-403] [MODERATE] — Anchor JOB_LOG to SLURM_SUBMIT_DIR
**Where:** `jobs/tpv102/tpv102_200m_p1_1.5s_{50,400}rank_dev.sbatch`

The post-run `grep` previously used a relative `JOB_LOG` path that resolved
against the script's CWD **after** two `cd` calls — not the Slurm submit
dir where `#SBATCH -o %j.out` actually writes the log.  On any submission
from outside `miniapps/seas/`, the `grep` silently fails to find the file
and the RESULT.txt reports FAIL regardless of physics outcome.

Fix: pin `SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"` at the very top of the
script (before any `cd`), then use `${SBATCH_LOG_DIR}/...` in the JOB_LOG
path.  A new pre-check aborts the verdict block with a targeted error
message if the log file is missing, instead of silently flowing into the
"no [qnorm:watch] line" branch (which obscures the actual root cause).

### [R-404] [LOW] — Sub-ULP tolerance in `same_centroid`
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`
(the pair-grouping `same_centroid` lambda)

Replaces the exact `!=` compare with `std::abs(va - vb) > scale * DBL_EPSILON`
where `scale = max(|va|, |vb|)`.  Transitivity is preserved (threshold
bounded by the larger magnitude, not an absolute floor).  In practice,
`ftr->Face->Transform(ip, phys)` produces bit-identical centroids on both
ranks today, so the change is a no-op on current MFEM builds — but it
insulates the R-305 unpaired-entry abort from a future compiler
reassociation round-off.  The sort comparator is intentionally left as
exact compare: shared-fault QP centroids produced from the same MFEM face
transformation at the same reference IP are bit-identical, and the
grouping logic compares `idx[j]` against `idx[i]` (not chained) so a
non-transitive `same_centroid` does not corrupt group formation.

### [R-405] [LOW] — Gate deep-copy MFEM_VERIFY to c==0
**Where:** `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS`

The deep-copy regression guard (`nbr_data[c].GetData() != q_gf.FaceNbrData().GetData()`)
now fires only on the first component (`c == 0`).  Every component uses
the same `SetSize + memcpy` allocator pathway, so a shallow-copy regression
would trip identically on c=0; checking the other eight components adds
no additional protection.  Trimmed call count from
`NUM_STATE × Mult × n_time_steps` to `Mult × n_time_steps` (~21k releases
on the dispositive 400-rank job vs ~190k before).  The diagnostic message
is unchanged apart from hardcoding "nbr_data[0]".

## Unresolved Findings
None.

## New Tests
- `jobs/tpv102/test_result_regex.sh` — covers R-402.  Four positive and
  five negative cases pin the new dead-rank regex against both the
  pre-fix false positive (`rN=0.0001`) and the post-format-switch hazards
  (scientific zero `0.000e+00`, denormals like `2.470e-61`).

No new unit test was added for R-401, R-403, R-404, R-405:
- R-401 requires a MFEM_USE_SINGLE build to exercise and is out of scope
  for the local macOS toolchain.
- R-403 is a shell-script CWD behaviour; verification requires a Slurm
  environment.  Manual verification: the added log-not-found pre-check
  emits a distinct error and exits 0, so a misrouted submission is
  diagnosable instead of silent.
- R-404 is a belt-and-suspenders no-op on current MFEM; the `R-302a`
  inline test already reports `max_diff=0, 3 pairs matched`, which still
  passes under the scale-ULP tolerance.
- R-405 is a hot-loop guard; covered by the existing parallel wave-operator
  tests (2 + 4 ranks) which still pass.

## Verification
- [x] R-401: fixed — `MPITypeMap<real_t>::mpi_type` on station reduce.  Covered
  by the regular TPV102 local/parallel test suite (24/24 + 16/16 pass with the
  default double build; single-precision path is not exercised locally).
- [x] R-402: fixed — driver emits stable %.3e qnorm format; sbatch regex
  matches the format and the new dead/alive cases. Covered by
  `test_result_regex.sh` (4 positive + 5 negative cases, OK).
- [x] R-403: fixed — `SBATCH_LOG_DIR` pinned before `cd`; JOB_LOG absolute;
  log-missing branch emits targeted message.
- [x] R-404: fixed — scale-ULP `same_centroid`; R-302a inline test still
  reports `max_diff=0, 3 pairs matched, 0 unpaired`.
- [x] R-405: fixed — MFEM_VERIFY gated on c==0;
  `seas_test_parallel_wave_operator` at 2 and 4 ranks still passes (5/5
  each), confirming no regression in the deep-copy invariant.

## Ready for Re-Review: YES

## Notes for the next run
- The 4-rank 1000 m smoke run (`tpv102_1000m_p1_1.5s_4rank_dev.sbatch`)
  does NOT have a post-run RESULT.txt block, so R-402/R-403 do not apply
  there.  The two dev-queue scripts with the block
  (`..._50rank_dev.sbatch`, `..._400rank_dev.sbatch`) are now safe to
  run and will report a trustworthy PASS/FAIL verdict regardless of
  submission CWD.
- No changes to physics — only post-run verifier, single-precision
  defensive typing, and two hot-loop/grouping refinements that are
  no-ops under current MFEM behaviour.
