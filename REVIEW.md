# Code Review: 2026-04-19 — v4 fresh review of TPV102 v3-fixes

This is a from-scratch adversarial review of the code after `tpv102_debug_v3_fix.md`
was applied. Each of the three passes was re-executed on the changed files; findings
below are not a re-verification of prior checklist items.

## Review Scope
- Plans / fix history consulted:
  - `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v1.md`
  - `tpv102_debug_v1_check.md` / `tpv102_debug_v1_fix.md`
  - `tpv102_debug_v2_check.md` / `tpv102_debug_v2_fix.md`
  - `tpv102_debug_v3_check.md` / `tpv102_debug_v3_fix.md`
- Files re-reviewed:
  - `miniapps/seas/drivers/tpv102_driver.cpp`
  - `miniapps/seas/dynamic/wave_operator.hpp`
  - `miniapps/seas/dynamic/wave_operator.inl`
  - `miniapps/seas/dynamic/tpv102_setup.hpp`
  - `miniapps/seas/io/paraview_output.hpp`
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp`
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch`
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch`
- Domain context:
  - `miniapps/seas/CLAUDE.md` (repo invariants, DG+fault physics rules)
  - MFEM `communication.hpp` (`MPITypeMap` specializations)
  - MFEM `pmesh.hpp` (shared-face / face-nbr / global-vertex APIs)
- Build + runtime verification:
  - `make seas_tpv102_driver` — **OK** on local macOS / conda mfem-dev
  - `mpirun -np 2 ./seas_test_r101_shared_fault` — **8/8 PASS**
    (R-101a / R-101b SKIP on TPV102 coarse mesh as disclosed; R-302a
    inline 2-tet test passes with `max_diff=0`, `3 pairs matched`)
  - `mpirun -np 4 ./seas_tpv102_driver --mesh tpv102_1000m.msh --tfinal 0.1 --debug-qnorm`
    — runs cleanly, no aborts, qnorm grows on all 4 ranks.

## Findings

### [R-401] [MODERATE] `dynamic/tpv102_setup.hpp:246-247` — TPV102StationWriter::Open (MPI variant) uses `MPI_DOUBLE` on a `std::vector<real_t>` buffer; silent memory corruption on `MFEM_USE_SINGLE` builds

**Category:** BUG (same root cause as R-303, missed in v3 fix scope)

**Description:**
R-303 fixed four `MPI_DOUBLE`-on-`real_t` sites: three in `tpv102_driver.cpp` and
one in `wave_operator.inl` (ctor h_min_). But `tpv102_setup.hpp` — which is
`#include`'d by `tpv102_driver.cpp` and compiled into `seas_tpv102_driver` —
still contains the same bug class on the station-ownership reduction:

```cpp
// dynamic/tpv102_setup.hpp, lines 230-247
std::vector<real_t> local_dist(nstations, std::numeric_limits<real_t>::max());
for (int s = 0; s < nstations; s++) { ... local_dist[s] = std::sqrt(...); }
std::vector<real_t> global_min_dist(nstations);
MPI_Allreduce(local_dist.data(), global_min_dist.data(), nstations,
              MPI_DOUBLE, MPI_MIN, comm);   // ← real_t buffer, MPI_DOUBLE
```

Both `local_dist` and `global_min_dist` are `std::vector<real_t>`. On a
`MFEM_USE_SINGLE` build (`real_t = float`, stride = 4 bytes), MPI_Allreduce
reads/writes 8 bytes per slot from/to a 4-byte-stride buffer. Result on single-
precision builds: garbage station-ownership distances → wrong ownership
tiebreaker → either duplicated or missing station files. On double-precision
builds (Frontera default, as far as the v3 fix doc implies) this is latent.

**Trigger:**
Any `MFEM_USE_SINGLE=YES` build that exercises TPV102. Does not block the
double-precision Frontera run but is a latent crash/corruption on the
supported single-precision config.

**Actual behavior:**
Works by coincidence on double builds (sizeof(real_t) == sizeof(double)).
Corrupts `global_min_dist[]` on single builds.

**Expected behavior:**
Use `MPITypeMap<real_t>::mpi_type`, matching the R-303 fix pattern already
applied in `tpv102_driver.cpp` and `wave_operator.inl`.

**Suggested fix:**
```diff
 // dynamic/tpv102_setup.hpp
-      // Global min distance across all ranks
-      std::vector<real_t> global_min_dist(nstations);
-      MPI_Allreduce(local_dist.data(), global_min_dist.data(), nstations,
-                    MPI_DOUBLE, MPI_MIN, comm);
+      // Global min distance across all ranks (R-401: match real_t type
+      // at compile time — same MPITypeMap pattern as R-303 fix).
+      std::vector<real_t> global_min_dist(nstations);
+      MPI_Allreduce(local_dist.data(), global_min_dist.data(), nstations,
+                    MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
```

**Test case:**
```python
@pytest.mark.skipif(not has_mfem_single_precision(),
                    reason="only applicable to MFEM_USE_SINGLE builds")
def test_R401_station_writer_single_precision():
    """
    On MFEM_USE_SINGLE build, TPV102StationWriter::Open must not
    produce corrupted station distances.  Precondition: real_t = float.
    Check no duplicated-station-file opens and all 9 stations get exactly
    one writer across the communicator.
    """
    run_tpv102_driver(nranks=4, build="single")
    station_writers = count_station_writers_across_ranks()
    # Exactly one rank opens each station file
    assert all(v == 1 for v in station_writers.values()), \
        f"duplicate/missing station owners under MFEM_USE_SINGLE: {station_writers}"
```

---

### [R-402] [MODERATE] `jobs/tpv102/tpv102_200m_p1_1.5s_{50,400}rank_dev.sbatch` — RESULT.txt regex false-positive marks alive-but-small ranks as FAIL; dispositive run's PASS/FAIL verdict can invert

**Category:** BUG (post-run verifier logic)

**Description:**
Both dev-queue sbatch scripts added in the R-306 fix use the same regex to
classify the final `[qnorm:watch]` line:

```bash
elif echo "${FINAL_WATCH}" | grep -qE "r[0-9]+=0(\b|[^.0-9])|r[0-9]+=1e-3[0-9]+"; then
   echo "FAIL: final [qnorm:watch] shows dead ranks..."
```

The first alternation `r[0-9]+=0(\b|[^.0-9])` is *not* a "dead rank" detector.
It matches any rank whose printed qnorm value starts with the digit `0`, followed
by anything non-digit-non-dot **or** a word boundary. In particular:

- `r10=0.0001` — rank 10 with ||Q||_∞ = 1e-4. MATCHES (`r10=0` + `\b` at the
  `0`→`.` transition since `0` is a word-char and `.` is not).
- `r10=0`    — true dead rank. MATCHES (correct behaviour).
- `r10=5e-21` — tiny but nonzero, doesn't start with `0`. Does **not** match.
- `r10=0.5`  — rank 10 with ||Q||_∞ = 0.5. MATCHES (false positive).

C++ `operator<<(ostream, double)` with default precision renders values in the
range `[1e-4, 1e6]` in *fixed* notation (e.g., `0.5`, `0.0001`), not in
scientific. In a working TPV102 run at tfinal=1.5 s, off-axis or edge ranks
can legitimately have ||Q||_∞ on the order of 0.1–1 depending on stress/velocity
scaling and output precision. Those ranks then trip the regex and the script
writes `FAIL: ...shows dead ranks...` even though the run was successful.

This is the verifier for the **dispositive 400-rank run**. A false FAIL means:
(a) node-hours are burned without a trustworthy verdict, (b) a follow-up human
eyeballs the log anyway, which was the problem R-306 was supposed to eliminate.

**Trigger:**
Any successful run whose final per-rank qnorm prints in fixed notation and
contains the leading digit `0`. Very likely to occur for edge ranks whose
wave perturbation is sub-unit but nonzero.

**Actual behavior:**
`r10=0.0001` → matches → FAIL written to `RESULT.txt`.
`r10=1e-21`  → no match → PASS written to `RESULT.txt` (correct).

**Expected behavior:**
The detector should recognize either (i) an exact literal zero (`=0` followed
by end-of-token), or (ii) a value in scientific notation with exponent ≤ -20
(the "essentially zero" floor). It must not match valid small fixed-notation
values.

**Suggested fix:**
Bound the detector to literal `=0` only at a token boundary (space or EOL),
and keep the existing scientific `1e-3XX` pattern:

```diff
-elif echo "${FINAL_WATCH}" | grep -qE "r[0-9]+=0(\b|[^.0-9])|r[0-9]+=1e-3[0-9]+"; then
+# R-402 fix: `r[0-9]+=0(\b|[^.0-9])` also matches valid small values like
+# `rN=0.0001`.  Restrict the "exact zero" case to `=0` followed by whitespace
+# or EOL, and keep the denormal-exponent case as-is.
+elif echo "${FINAL_WATCH}" | grep -qE "r[0-9]+=0([[:space:]]|$)|r[0-9]+=1e-3[0-9]+"; then
    echo "FAIL: final [qnorm:watch] shows dead ranks (rupture may not have crossed seam):" > "${RESULT_FILE}"
```

Belt-and-suspenders: also promote the driver to print qnorm with `%.3e`
formatting so the regex operates on a stable format instead of `operator<<`'s
value-dependent switching between fixed and scientific:

```diff
 // drivers/tpv102_driver.cpp, inside the [qnorm:watch] block
+std::cout << std::scientific << std::setprecision(3);
 for (int r : watch_unique)
 {
    std::cout << " r" << r << "=" << qn_all[r];
 }
 std::cout << " (hypo_rank=" << hypo_rank << ")\n";
+std::cout.unsetf(std::ios::scientific);    // restore default
```

**Test case:**
```bash
# jobs/tpv102/test_result_regex.sh
set -eu

# Positive cases: must match (FAIL)
for line in \
   "[qnorm:watch] r10=0 r11=1e5 (hypo_rank=10)"        \
   "[qnorm:watch] r10=1e-35 r11=1e5 (hypo_rank=10)"    \
   "[qnorm:watch] r10=1e-300 r11=1e5 (hypo_rank=10)"   \
; do
   echo "$line" | grep -qE "r[0-9]+=0([[:space:]]|$)|r[0-9]+=1e-3[0-9]+" \
      || { echo "REGRESS: should have matched: $line"; exit 1; }
done

# Negative cases: must NOT match (alive ranks)
for line in \
   "[qnorm:watch] r10=0.0001 r11=1e5 (hypo_rank=10)"   \
   "[qnorm:watch] r10=0.5 r11=1e5 (hypo_rank=10)"      \
   "[qnorm:watch] r10=5e-21 r11=1e5 (hypo_rank=10)"    \
   "[qnorm:watch] r10=1e6 r11=1e5 (hypo_rank=10)"      \
; do
   echo "$line" | grep -qE "r[0-9]+=0([[:space:]]|$)|r[0-9]+=1e-3[0-9]+" \
      && { echo "REGRESS: should NOT have matched: $line"; exit 1; }
   :
done

echo "R-402: RESULT.txt regex behaves correctly"
```

---

### [R-403] [MODERATE] `jobs/tpv102/tpv102_200m_p1_1.5s_{50,400}rank_dev.sbatch` — post-run `JOB_LOG` path is relative to a CWD the script itself changes; RESULT.txt silently records FAIL on any submission from outside `miniapps/seas/`

**Category:** BUG (CI/robustness)

**Description:**
Both scripts have the pattern (identical in 50-rank and 400-rank):

```bash
cd /scratch2/10024/zhaochun/seas-project/seas-mfem       # line ~53
# ...
cd miniapps/seas                                           # line ~60
make seas_tpv102_driver ...                                # line ~61
# ...
ibrun ./seas_tpv102_driver ...                             # line ~74
# ...
JOB_LOG="tpv102_200m_50r_${SLURM_JOB_ID}.out"              # relative name
FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)
```

Slurm writes `#SBATCH -o tpv102_200m_50r_%j.out` relative to the directory in
which `sbatch` was invoked (the submit CWD). The script then does
`cd ... && cd miniapps/seas`, so by the time the `grep "${JOB_LOG}"` runs, the
current working directory is `/scratch2/.../seas-mfem/miniapps/seas/`, which is
not the submit dir unless the user happened to invoke `sbatch` from exactly
`miniapps/seas/`. In every other submission workflow (e.g., submitting from
repo root, a login-node jump-host, or the `jobs/tpv102/` directory itself), the
`grep` returns "file not found" (exit 2), `FINAL_WATCH` ends up empty, and the
script then falls into the `elif [ -z "${FINAL_WATCH}" ]` branch and writes
`FAIL: no [qnorm:watch] line in job output — driver never produced diagnostic`
— even when the driver produced dozens of `[qnorm:watch]` lines successfully.

**Trigger:**
`sbatch jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch` from the repo root,
which matches the recommendation in the sbatch comments of rec #4 in the v3
check document:
> 4. `tpv102_200m_p1_1.5s_400rank_dev.sbatch` (~1.5 hr within dev-queue).
> **THE dispositive test.**

A user submitting from `/scratch2/10024/zhaochun/seas-project/seas-mfem/`
(the natural root to stand in for a repo root workflow) will see **every**
post-run RESULT.txt report FAIL regardless of physics outcome. R-306 fix's
stated goal ("PASS/FAIL machine-checkable verdict") becomes a no-op.

**Actual behavior:**
`grep "[qnorm:watch]" "tpv102_200m_400r_v2_${ID}.out" | tail -1` with CWD =
`miniapps/seas/` and log file at `/scratch2/.../seas-mfem/tpv102_200m_400r_v2_${ID}.out`
fails to find the file. FINAL_WATCH empty → FAIL verdict.

**Expected behavior:**
`JOB_LOG` must be an absolute path (or at minimum a path that resolves
independently of the script's own CWD changes).

**Suggested fix:**
Capture the submit CWD before any `cd`, and anchor the log path there:

```diff
 export LC_ALL=C
 export LANG=C

+# R-403 fix: pin SLURM_SUBMIT_DIR before we `cd` anywhere so the post-run
+# RESULT.txt check can locate the SBATCH -o file (Slurm writes it relative
+# to the submit dir, but the script changes CWD before grepping).
+SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
+
 module load intel/19.1.1
 # ...

 cd /scratch2/10024/zhaochun/seas-project/seas-mfem
 # ...
 cd miniapps/seas
 # ...
 ibrun ./seas_tpv102_driver ...
 # ...
 RESULT_FILE="${RESULT_DIR}/RESULT.txt"
-JOB_LOG="tpv102_200m_400r_v2_${SLURM_JOB_ID}.out"
+JOB_LOG="${SBATCH_LOG_DIR}/tpv102_200m_400r_v2_${SLURM_JOB_ID}.out"
+if [ ! -f "${JOB_LOG}" ]; then
+   echo "FAIL: post-run check cannot locate job log (looked at '${JOB_LOG}'). "\
+        "Hint: SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-<unset>}'." > "${RESULT_FILE}"
+   cat "${RESULT_FILE}"
+   exit 0
+fi
 FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)
```

Apply the identical change in the 50-rank variant.

**Test case:** N/A (script-level change). Manually verifiable by
`sbatch --hold` from the repo root, waiting for output file, then sourcing
the RESULT-writing block and asserting `grep FAIL "${RESULT_FILE}"` is empty
when the log contains a valid `[qnorm:watch]` line.

---

### [R-404] [LOW] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — centroid exact-equality pair matching silently skips the only-one-rank-hit-the-QP case on some topologies (corollary to R-305)

**Category:** POSSIBLE EDGE_CASE

**Description:**
R-305 (v3) promoted `n_unpaired > 0` to a hard abort. Good. But the centroid
pairing logic has one remaining soft edge: the `same_centroid` lambda uses
exact double equality on all 3 coordinates after `ftr->Face->Transform(ip, phys)`.
The v3 fix doc asserts (correctly) that both ranks compute the centroid from
the same MFEM face transformation at the same reference IP → bit-identical.

However, this bit-identity depends on an undocumented invariant: that MFEM's
`FaceElementTransformations::Face` returned by `GetSharedFaceTransformations`
on rank A and rank B, when driven with the same ip, produces bit-identical
`phys` vectors. If (a) the two ranks' face-nbr ghost vertex data differs in
any floating-point round-off on the face's geometric basis (it shouldn't;
Mesh::face_nbr_vertices is exchanged verbatim); or (b) an optimization pass in
the compiler reassociates float ops differently at the two call sites; then
two "same" QPs appear with different centroids and get classified as
**unpaired** on both sides, tripping R-305's abort on a false positive.

The current v3 test (R-302a inline) happens to have `max_diff=0` and
`3 pairs matched, 0 unpaired`, which argues the invariant holds on the local
mpich build. It is not an invariant MFEM documents, though.

**Suggested fix:** make `same_centroid` tolerant at the sub-ULP level only
(not the coordinate-scale tolerance that caused R-304). `1 ULP * |value|` is a
safe floor and preserves transitivity:

```diff
 auto same_centroid = [&](int a, int b)
 {
+   // Sub-ULP safety: `ftr->Face->Transform(ip, phys)` should produce bit-
+   // identical output on both ranks (MFEM exchanges face_nbr_vertices
+   // verbatim; same reference IP ⇒ same affine combination), but that is
+   // not a documented MFEM contract.  Allow <=1 ULP per coordinate against
+   // the coordinate's own magnitude so a future compiler-reassociation
+   // round-off doesn't flip matched pairs into R-305 aborts.
    for (int k = 0; k < 3; k++)
    {
-      if (all_data[a*REC + k] != all_data[b*REC + k]) { return false; }
+      double va = all_data[a*REC + k], vb = all_data[b*REC + k];
+      double scale = std::max(std::abs(va), std::abs(vb));
+      double ulp   = scale * std::numeric_limits<double>::epsilon();
+      if (std::abs(va - vb) > ulp) { return false; }
    }
    return true;
 };
```

Transitivity is preserved because the threshold is bounded by the *larger*
magnitude (not an absolute floor). This is a belt-and-suspenders change; the
current exact compare is defensible too.

**Test case:**
```python
def test_R404_same_centroid_ULP_safe():
    """1 ULP perturbation in one coord of one side must still match."""
    a = [1.0, 0.0, 0.0, 5.0]   # centroid + one field
    b = [np.nextafter(1.0, 2.0), 0.0, 0.0, 5.0]
    assert same_centroid(a, b) == True
```

---

### [R-405] [LOW] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — MFEM_VERIFY for nbr_data deep-copy (R-302 Part B) fires per-Mult call per component, adding O(NUM_STATE × Mult) release-build overhead for no useful protection

**Category:** QUALITY

**Description:**
```cpp
// wave_operator.inl:700-726
for (int c = 0; c < NUM_STATE; c++)
{
   for (int i = 0; i < ndof_total_; i++) { q_gf[i] = Q_data[...]; }
   q_gf.ExchangeFaceNbrData();
   const Vector &src = q_gf.FaceNbrData();
   nbr_data[c].SetSize(src.Size());
   std::memcpy(nbr_data[c].GetData(), src.GetData(), src.Size() * sizeof(real_t));
   MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(), ...);
}
```

R-302 Part B promoted this guard from `MFEM_ASSERT` to `MFEM_VERIFY` so the
invariant is enforced in Release. Intent understood — the comment calls the
deep-copy "load-bearing". But the assertion compares two pointers that are
architecturally guaranteed to differ immediately after `SetSize(n)` (which
allocates new storage via `Vector::SetSize`) followed by `std::memcpy`. The
check cannot fail *unless* a future refactor replaces `SetSize + memcpy` with
a shallow `operator=` — exactly the H2 scenario. In other words: it protects
against a refactor, not against a runtime condition the current code can
reach.

The check is cheap (two pointer reads + one compare) but runs
`NUM_STATE × n_Mult × n_time_steps` times on every ParMesh run of the driver;
~9 × 4 × 5250 = ~190k branches on the dispositive 400-rank job. Negligible
compute. QUALITY only — flagging it so the next reader understands the check's
intent is regression-catch, not runtime-safety.

**Suggested fix:** keep `MFEM_VERIFY` for the first iteration, downgrade
subsequent iterations to `MFEM_ASSERT` (debug-only):

```diff
+// R-302B reduced: one-time deep-copy invariant check on the first component;
+// subsequent components use the same allocator pathway so a single check
+// suffices at release.
 for (int c = 0; c < NUM_STATE; c++)
 {
    ...
    nbr_data[c].SetSize(src.Size());
    std::memcpy(nbr_data[c].GetData(), src.GetData(), src.Size() * sizeof(real_t));
-   MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
-               "nbr_data[" << c << "] aliases q_gf.FaceNbrData() — H2 ...");
+   if (c == 0)
+   {
+      MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
+                  "nbr_data[0] aliases q_gf.FaceNbrData() — H2 regression.  "
+                  "SetSize/memcpy path must be preserved in release.");
+   }
 }
```

**Test case:** N/A (quality refactor; coverage already via
`seas_test_parallel_wave_operator`).

---

## Summary
- Critical issues: **0**
- Moderate issues: **3**   (R-401 MPI datatype in station writer;
                            R-402 RESULT.txt regex false positive;
                            R-403 JOB_LOG path is CWD-dependent)
- Low issues: **2**        (R-404 sub-ULP same_centroid; R-405 hot-loop verify)
- Plan compliance: **FULL** — every R-301…R-309 item from
  `tpv102_debug_v3_check.md` is correctly and testably addressed, and the new
  `TestR302_InlineTwoTetSharedFault` exercises the shared-fault path end-to-end
  with `max_diff=0`.
- **Verdict: PASS WITH FIXES** — no finding in this review blocks the Frontera
  dispositive run from executing; the simulation code itself is sound and the
  R-005 / R-001 / R-101 fixes are in force. The two MODERATE sbatch issues are
  in the **post-run verifier** (RESULT.txt logic), not in the physics — but
  they will cause the dispositive run to be flagged FAIL even if it succeeds,
  which reintroduces exactly the "exit 0 ⇒ success" blind spot R-306 was
  supposed to close.
  Recommended action before burning 400-rank node-hours: apply R-402 + R-403
  fixes to both sbatch scripts so the post-run PASS/FAIL verdict is
  trustworthy. R-401 is advisory (only affects MFEM_USE_SINGLE builds, which
  the Frontera run does not use).

## Unreviewed Areas
- **Frontera module-load environment** — the sbatch scripts load
  `intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15
  fftw3/3.3.8`. This review does not re-verify the environment against a
  recent working Frontera build; if a module was retired or renamed between
  v2 and v3, the build step (`make seas_tpv102_driver`) fails before the
  driver ever runs. Mitigation: dry-run the module loads on a Frontera login
  node before submitting dev-queue jobs.
- **METIS partitioning of the 200 m mesh onto 400 ranks** — whether the 400-
  rank partition produces a fault-face-owning rank whose neighbors see bulk-
  wave energy arriving within ~0.5 s depends on METIS's placement, which is
  undocumented. This review assumes it behaves as the 4-rank smoke test
  observed.
- **`MFEM_USE_SINGLE` build** — build was not attempted in single precision;
  R-401 is inferred from code shape, not observed.
- **R-101 runtime diagnostic on a mesh with non-duplicated-vertex fault** —
  R-302a inline test (2 tets) passes. No larger non-duplicated-vertex mesh
  exists in the repo, so the diagnostic's behaviour under realistic QP counts
  (e.g., 10^4 shared fault QPs) is unexercised.
