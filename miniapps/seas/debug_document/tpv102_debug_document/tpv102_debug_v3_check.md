# Code Review v3: Fresh adversarial review of v2 fix round

Re-auditing the code *after* `tpv102_debug_v2_fix.md` was applied.  This
is a from-scratch review — every pass re-executed.

## Review Scope
- Input docs:
  - `tpv102_debug_v1.md`
  - `tpv102_debug_v1_check.md` / `tpv102_debug_v1_fix.md`
  - `tpv102_debug_v2_check.md` / `tpv102_debug_v2_fix.md`
- Files reviewed (diffed against `git HEAD`):
  - `miniapps/seas/drivers/tpv102_driver.cpp`
  - `miniapps/seas/dynamic/tpv102_setup.hpp`
  - `miniapps/seas/dynamic/wave_operator.hpp`
  - `miniapps/seas/dynamic/wave_operator.inl`
  - `miniapps/seas/io/paraview_output.hpp`
  - `miniapps/seas/Makefile`
- Added files reviewed:
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp`
  - `miniapps/seas/jobs/tpv102/tpv102_1000m_p1_1.5s_4rank_dev.sbatch`
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch`
  - `miniapps/seas/jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch`
  - `miniapps/seas/jobs/tpv102/tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch`
- Domain context:
  - `miniapps/seas/CLAUDE.md` (project invariants)
  - `miniapps/seas/dynamic/fault_face_flux.cpp` (Evaluate signature and fields mutated)
  - MFEM `pmesh.hpp` (public face-neighbour API)

## Findings

### [R-301] [CRITICAL] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — Only 3 of the 8 mutable `DOFData` fields are checked; divergence in `V2`, `tau2_corr`, `sigma_n_corr`, `slip1`, `slip2` is silently accepted

**Category:** BUG

**Description:**
`VerifySharedFaultDOFDataConsistency` is the runtime guard that replaces the "unverified MFEM invariant" from the v2 review. Its whole purpose is to assert bit-equality of `DOFData` between the two ranks that own every shared fault QP. The implementation packs a 6-double record per QP (lines 1054-1059):

```cpp
local_data.push_back(phys(0));         // centroid x
local_data.push_back(phys(1));         // centroid y
local_data.push_back(phys(2));         // centroid z
local_data.push_back(d.tau1_corr);     // ← checked
local_data.push_back(d.V1);            // ← checked
local_data.push_back(d.psi);           // ← checked
```

`DOFData` mutable state that `FaultFaceFlux::Evaluate` writes is **eight** fields (see `fault_face_flux.cpp:150-155`):

```cpp
data.slip_rate = V_abs;
data.V1 = V1;
data.V2 = V2;                      // ← NOT CHECKED
data.tau1_corr = ...;
data.tau2_corr = ...;              // ← NOT CHECKED
data.sigma_n_corr = ...;           // ← NOT CHECKED
// plus driver-level RK4 averaging writes slip1, slip2, slip_rate
data.slip1 += dof_data[i].V1 * dt_step;   // ← NOT CHECKED
data.slip2 += dof_data[i].V2 * dt_step;   // ← NOT CHECKED
data.slip_rate = std::sqrt(V1² + V2²);    // ← NOT CHECKED
```

TPV102 is a nominally mode-II fault (strike-slip), so `V2`/`tau2_corr` are small compared to `V1`/`tau1_corr` — but "small" is exactly the regime where a subtle orientation sign-flip is hardest to see in `V1` and easiest to see in `V2`. In a mixed-mode setup (TPV101, TPV102 with off-axis nucleation, or any non-planar fault), this omission silently hides the very bug the diagnostic is supposed to catch.

Combined with R-302 below, this means the v2 fix ships a "runtime guard" that is not an actual guard for the full invariant it is documented to enforce.

**Trigger:**
Any build where R-001's (+,−) swap is active (i.e., the mesh produces shared fault faces). The bug is latent but flagged at MODERATE only because current TPV102 meshes bypass the shared-fault path (see R-302).

**Actual behavior:**
```cpp
constexpr int REC = 6;                         // ← should be 11
local_data.push_back(phys(0));
…
local_data.push_back(d.psi);                   // stops here
```

`max_diff_field == 3|4|5` mapped to `"tau1_corr"|"V1"|"psi"`. Anything else is silently labelled "psi" in the abort message (default-case fallthrough — see R-307).

**Expected behavior:**
Pack every field that `Evaluate` or the driver's RK4-averaging step mutates. At minimum: `tau1_corr, tau2_corr, sigma_n_corr, V1, V2, psi, slip1, slip2`.

**Suggested fix:**
```diff
-      // Pack per-QP records: (cx, cy, cz, tau1_corr, V1, psi).
-      // One record = 6 doubles.  Each rank contributes one record per
-      // shared fault QP it owns; the same physical QP appears on both
-      // ranks that share the fault face.
-      constexpr int REC = 6;
+      // Pack per-QP records: centroid + every DOFData field that
+      // FaultFaceFlux::Evaluate or the driver-side RK4 averaging writes.
+      // One record = 11 doubles.
+      constexpr int REC = 11;
+      constexpr int FIELD_BASE = 3;   // index in record where fields start
+      static const char* FIELD_NAMES[8] = {
+         "tau1_corr", "tau2_corr", "sigma_n_corr",
+         "V1", "V2", "psi", "slip1", "slip2"
+      };
…
            local_data.push_back(phys(0));
            local_data.push_back(phys(1));
            local_data.push_back(phys(2));
            local_data.push_back(d.tau1_corr);
+           local_data.push_back(d.tau2_corr);
+           local_data.push_back(d.sigma_n_corr);
            local_data.push_back(d.V1);
+           local_data.push_back(d.V2);
            local_data.push_back(d.psi);
+           local_data.push_back(d.slip1);
+           local_data.push_back(d.slip2);
…
-      const char *field_name = (max_diff_field == 3) ? "tau1_corr"
-                              : (max_diff_field == 4) ? "V1"
-                              : "psi";
+      const char *field_name = (max_diff_field >= FIELD_BASE &&
+                                max_diff_field <  FIELD_BASE + 8)
+                             ? FIELD_NAMES[max_diff_field - FIELD_BASE]
+                             : "<unknown>";
```

**Test case:**
```python
def test_R301_verifier_checks_all_mutable_fields():
    """
    Stub DOFData so that V2 drifts by 1 ULP between ranks while V1,
    tau1_corr, psi stay identical.  Expect MFEM_ABORT; current code
    would silently accept.
    """
    with patch_evaluate_inject_v2_drift():
        with pytest.raises(RuntimeError, match="V2"):
            run_r101_test(nranks=2, mesh="two_tet_shared_fault.msh")
```

---

### [R-302] [CRITICAL] `tests/parallel/test_r101_shared_fault.cpp` — SKIPS on every current TPV102 mesh; R-001/R-101/R-102 are shipped as defensive code for a path that does not run in production, with zero runtime coverage

**Category:** DEVIATION / ASSUMPTION

**Description:**
`tpv102_debug_v2_fix.md` § "R-101 unit-test caveat" acknowledges:

> "MFEM does NOT classify TPV102-style fault faces (generated by Gmsh's
>  `BooleanFragments`, which duplicates vertices on either side of the
>  fault surface) as *shared* faces after `ParMesh` partitioning.  Each
>  rank's portion of the fault appears in its `GetNSharedFaces()` count
>  as zero and in its `GetNBE()` list as local interior-BE faces."
>
> "The `ComputeSharedFaceFluxRHS` shared-fault code path is NOT reached
>  for these meshes — the fault is handled entirely by
>  `ComputeFaceFluxRHS`'s interior-fault branch."
>
> "R-001's (+,-) swap is defensive code for a code path that never
>  fires on the current meshes, and thus cannot be the cause of the v1
>  'rupture-stops-at-partition-seam' symptom."

This is a factual disclosure in the fix report — but it also means the test and the diagnostic are *vacuous* on every mesh in the repo, and the `tpv102_debug_v2_fix.md` "Test suite: PASS — `seas_test_r101_shared_fault` (2 ranks) — SKIPPED (see below)" entry is a *false* pass signal. SKIPPED ≠ PASS when the whole point of the test is to verify the fix.

Concretely:
- **No runtime path in current production** exercises `ComputeSharedFaceFluxRHS`'s shared-fault branch.
- **No test** ever activates `VerifySharedFaultDOFDataConsistency`'s body — the `any_shared_global == 0` short-circuit at line 175-180 returns early.
- The **actual** fix for the v1 H1 symptom is almost certainly R-005 (the explicit memcpy that replaced the aliased `nbr_data[c] = q_gf.FaceNbrData()` assignment). But R-005 is guarded only by `MFEM_ASSERT` (debug-only), so in Release builds the pointer-alias invariant is silently trusted. If a future refactor replaces the `SetSize + memcpy` with any form of `operator=` that happens to take a shallow branch, the v1 bug returns undetected.

**Trigger:**
Anyone who reads the check doc or the fix report and concludes "the v1 bug has been fixed and regression-tested" is acting on false assurance.

**Actual behavior:**
- `ComputeSharedFaceFluxRHS` shared-fault block contains R-001 swap + R-102 guards; none are reachable on any mesh in the repo.
- `VerifySharedFaultDOFDataConsistency` returns early on every TPV102 partitioning tried so far.
- The 2-rank MPI unit test runs both `TestR101_OneStage` and `TestR101_TenSteps` and reports "SKIPPED" for both. The Makefile target `test-r101-shared-fault` completes with exit 0, which aggregators read as "PASS".

**Expected behavior:**
- Either (a) augment the mesh so the shared-fault path is exercisable, so the test actually verifies the fix; or (b) explicitly mark the test's skip as INCONCLUSIVE (non-zero exit, or at minimum a CI-visible warning), so the claim "R-001 is verified" cannot be made from the automated suite.
- The actual fix (R-005) should be promoted from `MFEM_ASSERT` to `MFEM_VERIFY` so the alias invariant is checked in Release, since that is the load-bearing fix.

**Suggested fix (two-part):**

Part A — exercise the shared-fault path:

```diff
 # tests/parallel/test_r101_shared_fault.cpp
+// Build a 2-tet mesh with a shared fault face that MFEM actually
+// classifies as "shared" (not duplicated-vertex).  Requires an inline
+// mesh constructor that doesn't use BooleanFragments:
+static std::string BuildTwoTetSharedFaultMeshInline(int rank);
```

Shipping a 2-tet inline mesh where the shared face is a true shared face is the only way to make this test meaningful. A simple way: build a 2-tet mesh programmatically with `Mesh::Mesh(dim, nv, nelem, nbdry)` rather than reading a Gmsh file; partition so each tet goes to a different rank; the inter-tet face will be a real shared face with `ParMesh::GetNSharedFaces() == 1`.

Part B — promote R-005 from `MFEM_ASSERT` to `MFEM_VERIFY`:

```diff
 const Vector &src = q_gf.FaceNbrData();
 nbr_data[c].SetSize(src.Size());
 std::memcpy(nbr_data[c].GetData(), src.GetData(),
             src.Size() * sizeof(real_t));
-MFEM_ASSERT(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
-            "nbr_data[c] aliases q_gf.FaceNbrData() — H2 regression");
+MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
+            "nbr_data[" << c << "] aliases q_gf.FaceNbrData() — H2 "
+            "regression.  The deep-copy invariant is load-bearing for "
+            "cross-rank bulk wave propagation and must hold in Release.");
```

**Test case:**
```python
def test_R302_inline_two_tet_shared_fault_actually_exercises_swap():
    """Build a 2-tet inline mesh with a real MFEM shared face tagged as
    fault.  Expect GetNSharedFaults() == 1 on each rank, and expect
    VerifySharedFaultDOFDataConsistency to NOT short-circuit."""
    ctx = build_two_tet_shared_fault_context(nranks=2)
    assert ctx.num_shared_fault > 0, "test must hit the shared-fault path"
    assert not ctx.verify_returned_early(), \
        "VerifySharedFaultDOFDataConsistency must actually run its body"
```

---

### [R-303] [CRITICAL] `drivers/tpv102_driver.cpp` — `--debug-qnorm` gather uses `MPI_DOUBLE` on a `std::vector<real_t>` buffer; memory corruption on `MFEM_USE_SINGLE` builds

**Category:** BUG

**Description:**
```cpp
real_t qn_local = Q.Normlinf();
std::vector<real_t> qn_all;
if (rank == 0) { qn_all.resize(nprocs); }
MPI_Gather(&qn_local, 1, MPI_DOUBLE,
           rank == 0 ? qn_all.data() : nullptr, 1, MPI_DOUBLE,
           0, comm);
```

`qn_local` is `real_t`.  `qn_all` stores `real_t`.  But the MPI datatype is hardcoded `MPI_DOUBLE`.  With the standard `real_t = double` build this is coincidentally correct; with `MFEM_USE_SINGLE` enabled (`real_t = float`, `sizeof(real_t) = 4`), MPI_Gather reads/writes 8 bytes per slot from/to a 4-byte-stride buffer.  This is **silent memory corruption**: ranks past rank 0 work fine, rank 0 writes `qn_all[r] * 2` slots with stale-aligned data.

The MFEM-idiomatic fix is to use `MPITypeMap<real_t>::mpi_type` (MFEM's template) or the compile-time-selected `MPI_REAL_T`.

**Trigger:**
Any build with `MFEM_USE_SINGLE = YES`. Also any code-review tool that notices the type mismatch.

**Actual behavior:**
Works silently on double builds; corrupts `qn_all` on single-precision builds.

**Expected behavior:**
MPI datatype must match buffer element type.

**Suggested fix:**
```diff
 real_t qn_local = Q.Normlinf();
 std::vector<real_t> qn_all;
 if (rank == 0) { qn_all.resize(nprocs); }
-MPI_Gather(&qn_local, 1, MPI_DOUBLE,
-           rank == 0 ? qn_all.data() : nullptr, 1, MPI_DOUBLE,
-           0, comm);
+const MPI_Datatype mpi_real_t =
+   (sizeof(real_t) == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT;
+MPI_Gather(&qn_local, 1, mpi_real_t,
+           rank == 0 ? qn_all.data() : nullptr, 1, mpi_real_t,
+           0, comm);
```

(Preferred: `MPITypeMap<real_t>::mpi_type` if accessible here; the sizeof-based conditional is a minimal-deps fallback.)

**Test case:**
```python
@pytest.mark.skipif(not has_mfem_single_precision(),
                    reason="only applicable to MFEM_USE_SINGLE builds")
def test_R303_qnorm_gather_correct_dtype():
    """
    On a single-precision MFEM build, the --debug-qnorm output must not
    show memory-corruption patterns (NaN, 1e-300 garbage) on non-root
    ranks.  Precondition: MFEM_USE_SINGLE=YES; real_t = float.
    """
    log = run_driver_with_debug_qnorm(nranks=4)
    qnorm_values = parse_debug_qnorm(log)
    for v in qnorm_values:
        assert not math.isnan(v)
        assert abs(v) < 1e100      # detect garbage-padded doubles
```

---

### [R-304] [MODERATE] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — `std::sort` comparator is **not transitive** across coordinate tolerance; formal UB and produces inconsistent grouping near-tied centroids

**Category:** BUG (undefined behavior)

**Description:**
```cpp
std::sort(idx.begin(), idx.end(), [&](int a, int b) {
   for (int k = 0; k < 3; k++) {
      double va = all_data[a*REC + k], vb = all_data[b*REC + k];
      if (std::abs(va - vb) > coord_tol) { return va < vb; }
   }
   return false;
});
```

`std::sort` requires **strict weak ordering**: if `!(a<b) && !(b<a)` then `a` and `b` are equivalent and by transitivity any `c` equivalent to `a` is equivalent to `b`. This comparator violates that:

- `A = (0, 0, 0)`, `B = (0.5e-6, 0, 0)`, `C = (1.5e-6, 0, 0)`, with `coord_tol = 1e-6`.
- `A vs B`: |0 − 0.5e-6| = 0.5e-6 ≤ tol → neither < the other.  A ≈ B.
- `B vs C`: |0.5e-6 − 1.5e-6| = 1e-6 ≤ tol (strict `>`) → neither < the other.  B ≈ C.
- `A vs C`: |0 − 1.5e-6| = 1.5e-6 > tol → A < C.

So A ≈ B ≈ C but A < C — violates transitivity. Behavior of `std::sort` is **undefined** under violations of strict weak ordering; on some libc++ implementations this can crash (heap corruption).

In practice TPV102 QP centroids are separated by at least the face edge length (~200 m at 200 m mesh), so the tolerance ambiguity does not trigger. But as soon as someone passes a finer tolerance or refines to h < 1e-6 m (impossible for us today, but the code doesn't know that), this becomes a latent crash.

**Trigger:**
Any mesh where two distinct shared-fault QP centroids are within `coord_tol` of each other on one axis but not on all axes. Currently impossible on any real mesh, but a latent bug waiting for a tolerance tweak.

**Actual behavior:**
UB.  Most STLs happen to handle it gracefully, but that is not guaranteed.

**Expected behavior:**
Either (a) use an exact lexicographic comparator (no tolerance, fine since both ranks produce identical coordinates via `ftr->Face->Transform(ip)` with the same reference IP on the same MFEM Face), or (b) quantize coordinates to a grid and compare exactly.

**Suggested fix:**
Replace the tolerance-based comparator with lexicographic exact compare.  Shared fault QPs on the two ranks that own them **are** computed from the same face transformation evaluated at the same reference IP → bit-identical.  The tolerance was unnecessary.

```diff
 std::sort(idx.begin(), idx.end(), [&](int a, int b)
 {
-   for (int k = 0; k < 3; k++)
-   {
-      double va = all_data[a*REC + k], vb = all_data[b*REC + k];
-      if (std::abs(va - vb) > coord_tol) { return va < vb; }
-   }
-   return false;
+   // Exact lex compare: both ranks produce bit-identical centroids
+   // from the same mesh face geometry + same reference IP.
+   for (int k = 0; k < 3; k++)
+   {
+      double va = all_data[a*REC + k], vb = all_data[b*REC + k];
+      if (va != vb) { return va < vb; }
+   }
+   return false;
 });
…
 auto same_centroid = [&](int a, int b)
 {
    for (int k = 0; k < 3; k++)
    {
-      if (std::abs(all_data[a*REC + k] - all_data[b*REC + k]) > coord_tol)
-      { return false; }
+      if (all_data[a*REC + k] != all_data[b*REC + k]) { return false; }
    }
    return true;
 };
```

**Test case:**
```python
def test_R304_sort_comparator_strict_weak_ordering():
    """
    Synthetic test: 3 colinear centroids that straddle the tolerance.
    Current comparator gives UB; fix makes sort stable.
    """
    # A, B, C along x: A=0, B=coord_tol/2, C=coord_tol*1.5
    centroids = np.array([[0, 0, 0],
                          [0.5e-6, 0, 0],
                          [1.5e-6, 0, 0]])
    order = sort_with_current_comparator(centroids, tol=1e-6)
    # Strict weak ordering requires unique sort order; many STLs return
    # inconsistent ordering here.  Post-fix: deterministic lex order.
    order_fix = sort_with_exact_comparator(centroids)
    assert order_fix == [0, 1, 2]
```

---

### [R-305] [MODERATE] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — `group_size != 1 && != 2` accumulates `n_unpaired` but the function never reports or aborts on this anomaly; a non-2 group silently indicates a mesh-topology bug

**Category:** EDGE_CASE

**Description:**
```cpp
if (group_size == 2) { … n_pairs++; }
else if (group_size != 1) {
   // A shared face is owned by exactly 2 ranks; anything else means
   // the centroid-match collapsed unrelated faces, which would
   // already be a logic error in this diagnostic.  Still, record
   // for the diagnostic message.
   n_unpaired += group_size;
}
else { n_unpaired++; }
```

`n_unpaired` is computed, printed on success ("N pairs matched, n_unpaired entries, max_diff=… OK"), but the function's success condition (`fail_global == 0`) only checks `max_diff > tol` — it does not abort on `n_unpaired > 0`. The comment admits "would already be a logic error in this diagnostic," but the code does not abort; it prints an OK message with a nonzero `n_unpaired` count buried in it.

Concretely: if a shared fault QP shows up on **only one** rank (mesh topology bug, or an R-102-style peer-resolution failure that MFEM_VERIFY did not catch because `resolved` was true but `peer_rank` was stale), the centroid group has size 1, `n_unpaired++`, but the diagnostic reports success.

**Trigger:**
Any mesh-partitioning pathology that produces a fault QP owned by only one rank (or more than two).

**Actual behavior:**
Reports "OK" with nonzero `n_unpaired` in the message.

**Expected behavior:**
Treat `n_unpaired > 0` as a diagnostic failure.

**Suggested fix:**
```diff
 int fail_local = (max_diff > tol) ? 1 : 0;
+if (n_unpaired > 0) { fail_local = 1; }
 int fail_global = 0;
 MPI_Allreduce(&fail_local, &fail_global, 1, MPI_INT, MPI_MAX, comm);

 if (fail_global)
 {
+   if (n_unpaired > 0 && max_diff <= tol)
+   {
+      MFEM_ABORT("R-101 shared-fault DOFData: " << n_unpaired
+         << " unpaired entries (every shared QP should have exactly "
+         << "2 ranks).  Likely mesh-partitioning pathology or "
+         << "peer-rank resolution failure.");
+   }
    …
 }
```

**Test case:**
```python
def test_R305_unpaired_aborts():
    """If a shared fault QP is owned by only one rank, diagnostic must
    abort, not print 'OK'."""
    # Force unpair by post-processing dof_data on rank 1 to move its
    # shared QP's centroid before VerifySharedFaultDOFDataConsistency.
    with force_unpaired_shared_fault_qp():
        with pytest.raises(RuntimeError, match="unpaired entries"):
            run_driver_and_verify_r101()
```

---

### [R-306] [MODERATE] `jobs/tpv102/tpv102_200m_p1_1.5s_400rank_dev.sbatch` — No post-run machine-checkable success criterion; the dispositive run reports exit 0 regardless of whether the v1 symptom recurred

**Category:** QUALITY (CI)

**Description:**
The 400-rank dev-queue script is described in the fix report as "the DISPOSITIVE test" for the v1 fix. If it passes (rupture crosses the first partition seam), we proceed to production; if it fails, we do not. But the script's exit status is determined by `ibrun ./seas_tpv102_driver …` — which returns 0 on normal completion regardless of whether the physics outcome was correct. The expected signature ("rupture past |x| ≥ 1.5 km by t=1.5 s", or "`[qnorm:watch]` shows non-hypocenter ranks > 0") is documented in the sbatch comments but **not** checked programmatically.

A human must eyeball the job output and/or open ParaView. For a CI-driven workflow — or for "did my node-hours burn productively?" at 400 ranks — this is insufficient. The fix report's next-step workflow instructs "look for `[qnorm:watch]` showing …" but does not include a grep-based PASS/FAIL write-out.

**Trigger:**
Any human failure to read the output carefully.  Also: re-runs where the user assumes `exit 0` = success.

**Actual behavior:**
`ibrun` exits 0 on driver completion even if rupture stalled at the hypocenter.

**Expected behavior:**
After the driver completes, the sbatch script should parse the `[qnorm:watch]` final line and verify non-hypocenter ranks have non-trivial ||Q||_∞; parse station output for `V1 > V_ini`; write `PASS` / `FAIL` to a result file the user can immediately `grep`.

**Suggested fix:**
Add post-run validation to the 400-rank sbatch (and the 50-rank one):

```diff
 ibrun ./seas_tpv102_driver \
       --mesh tpv102/mesh/tpv102_200m.msh \
       …
       --debug-qnorm
+
+# Post-run PASS/FAIL check: non-hypocenter ranks must register
+# ||Q||_inf > 1e-10 by t=1.5 s.
+RESULT_FILE="${RESULT_DIR}/RESULT.txt"
+JOB_LOG="tpv102_200m_400r_v2_${SLURM_JOB_ID}.out"
+FINAL_WATCH=$(grep "\[qnorm:watch\]" "${JOB_LOG}" | tail -1)
+if echo "$FINAL_WATCH" | grep -qE "r[0-9]+=0\b|r[0-9]+=1e-3[0-9]+" ; then
+   echo "FAIL: final qnorm watch line shows dead ranks: $FINAL_WATCH" > "${RESULT_FILE}"
+elif grep -q "R-101 shared-fault DOFData consistency FAILED" "${JOB_LOG}"; then
+   echo "FAIL: R-101 diagnostic tripped"    > "${RESULT_FILE}"
+else
+   echo "PASS: cross-rank wave propagation evidence present"  > "${RESULT_FILE}"
+fi
+cat "${RESULT_FILE}"
```

**Test case:**  N/A — script-level change; verified by running it once on dev queue.

---

### [R-307] [LOW] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — Abort-message field label falls through to "psi" for any unknown index, muddying the diagnostic

**Category:** QUALITY

**Description:**
```cpp
const char *field_name = (max_diff_field == 3) ? "tau1_corr"
                        : (max_diff_field == 4) ? "V1"
                        : "psi";
```

If `max_diff_field == -1` (the sentinel value when no pair is found — admittedly unreachable given the current structure, but fragile to future edits) or any other unexpected index, the abort message mis-labels it as `psi`. Combined with R-301's fix for more fields, the cascade becomes unwieldy; use a table.

**Suggested fix:** See R-301 suggested fix (adds `FIELD_NAMES[]` table with a bounds check and `"<unknown>"` fallback).

**Test case:** covered by R-301's test.

---

### [R-308] [LOW] `io/paraview_output.hpp:CommitSchedule` — Takes `cycle` and `V_max` parameters that are unused (`/*cycle*/`, `/*V_max*/`); the trio `PeekShouldWrite → CommitSchedule` is easy to call with the wrong argument triple

**Category:** QUALITY

**Description:**
```cpp
void CommitSchedule(int /*cycle*/, real_t time, real_t /*V_max*/)
{
   last_write_time_ = time;
}
```

The signature mirrors `Save(cycle, time, V_max)` and `PeekShouldWrite(cycle, time, V_max)` for call-site symmetry, but `cycle` and `V_max` are discarded.  A reader must cross-reference `Save` to understand why three args exist when only one is used.  More concerning: a future caller might read this as `CommitSchedule(step, last_write_time_, 0)` and expect the last two args to matter.

**Suggested fix:**
```diff
-void CommitSchedule(int /*cycle*/, real_t time, real_t /*V_max*/)
-{
-   last_write_time_ = time;
-}
+/// Advance the internal write-time watermark (no cycle/V_max arg — the
+/// single-time-advance is the only thing CommitSchedule does).  Paired
+/// with `PeekShouldWrite(cycle, time, V_max)`; the caller is expected
+/// to invoke `PeekShouldWrite` first and `CommitSchedule` only once the
+/// per-cycle writes have completed.
+void CommitSchedule(real_t time) { last_write_time_ = time; }
```

Then update the driver:
```diff
-pv_out->CommitSchedule(step_num, time, V_max);
+pv_out->CommitSchedule(time);
```

**Test case:** N/A — refactor-only.

---

### [R-309] [LOW] `tests/parallel/test_r101_shared_fault.cpp:BuildFaultCoords` — No null check on `GetInteriorFaceTransformations` / `GetSharedFaceTransformations`; crashes if MFEM returns null

**Category:** EDGE_CASE

**Description:**
```cpp
for (int i = 0; i < int_faces.Size(); i++)
{
   auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[i]);
   const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
   …
}
```

`wave.GetFaultInteriorFaces()` only contains faces that passed a non-null check in the WaveOperator ctor, so `ftr` should be non-null in practice. But the test adds a second call to `GetInteriorFaceTransformations` outside the ctor; if MFEM's internal face bookkeeping ever becomes stateful (e.g., after a mesh refine), the cached list might have a stale face that now returns null. The test would crash with a segfault, not a clean FAIL.

**Suggested fix:**
```diff
 for (int i = 0; i < int_faces.Size(); i++)
 {
    auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[i]);
+   if (!ftr)
+   {
+      MFEM_ABORT("test: interior fault face " << int_faces[i]
+                 << " returned null FTR after WaveOperator ctor filtered it.  "
+                 "MFEM invariant violation — check mesh state.");
+   }
    const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
    …
 }
```

(Same pattern for the shared-faces loop below.)

**Test case:** covered by R-302's new inline-mesh test.

---

## Recommended additional unit tests (summary)

1. **Inline-mesh R-101 test** (R-302 Part A).  Build a 2-tet mesh programmatically with a real MFEM shared face; partition 1 tet per rank; confirm `GetNSharedFaces() == 1` and the shared-fault path actually executes.  This is the **only** way to verify R-001 without waiting for a mesh format that classifies fault faces as shared.
2. **Full-field DOFData consistency test** (R-301).  Inject a 1-ULP drift in `V2` (or `sigma_n_corr`, `tau2_corr`, `slip2`) via a patched Evaluate stub; confirm the diagnostic aborts.  Current diagnostic would silently accept.
3. **Comparator strict-weak-ordering test** (R-304).  Synthetic centroids that straddle the tolerance; expose the UB.  Validates the exact-compare fix.
4. **Unpaired-entry test** (R-305).  Remove one side of a shared fault QP (mock the mesh); confirm diagnostic aborts.
5. **`MFEM_USE_SINGLE` build sanity** (R-303).  Build with `MFEM_USE_SINGLE=YES`; run `--debug-qnorm`; assert no NaN/garbage in the qnorm output.  Lightweight CI pre-check.

## Recommended short Frontera dev runs (unchanged from v2, with priority ordering)

1. `tpv102_1000m_p1_1.5s_4rank_dev.sbatch` (~15 min).  Cheapest sanity; confirms cross-rank bulk-wave energy flows.  **Must pass before submitting 2 or 3.**
2. `tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch` (~15 min, optional visual smoke).  Open the fault-surface PVD in ParaView; look for partition-seam stripes in `slip_rate_strike`.  If stripes visible → R-001 fix incomplete (or, given R-302, R-005 fix incomplete).
3. `tpv102_200m_p1_1.5s_50rank_dev.sbatch` (~2 hr).  Intermediate scale; exercise R-005 deep-copy at ~50 graded-partition seams.  **Add the R-306 post-run PASS/FAIL write-out before submitting.**
4. `tpv102_200m_p1_1.5s_400rank_dev.sbatch` (~1.5 hr within dev-queue).  **THE dispositive test.**  Do not burn node-hours here until 1+2+3 pass.  **Requires R-306 post-run check** so a regressed run doesn't hide behind `exit 0`.

### New recommendation v3: a zero-knowledge version of dev-run #4

Before committing to the 400-rank run, add a 50-rank variant that runs for exactly the nucleation window (0.0–1.2 s) and greps for the *positive* signal ("`[qnorm:watch]`" shows non-hypocenter-rank ||Q||_∞ growing monotonically). This is structurally the same as rec #3 but halved in physical time, so it fits in ~1 hr and pre-flights the dispositive run.  If it fails, the 400-rank run is pointless; if it passes, the 400-rank run is high-confidence.

## Summary
- Critical issues: 3 (R-301 partial-field diagnostic, R-302 untested defensive code for dead path, R-303 MPI datatype mismatch on single-precision)
- Moderate issues: 3 (R-304 non-transitive sort comparator, R-305 unpaired accumulator doesn't abort, R-306 no machine-checkable success criterion on dispositive run)
- Low issues: 3 (R-307 abort-message label cascade, R-308 CommitSchedule unused args, R-309 test null-check)
- Plan compliance: **PARTIAL**.  The v2 fix report ticks every v2 finding, but two of them (R-101, R-102) are defensive code for a path that never fires on the real meshes; the actual fix (R-005) is guarded only by a debug-mode assert; the unit test SKIPs on every real mesh.
- Verdict: **PASS WITH FIXES.**  R-301 (diagnostic missing fields) and R-302 (test vacuous) are the blockers for trusting the claim that R-001/R-101/R-102 "work" — they may, but we have no evidence.  R-303 is a latent crash on a supported build config.  R-304/R-305/R-306 should be fixed before the 400-rank run so its outcome is actually interpretable.

## Unreviewed Areas
- Actual R-005 efficacy under Release builds: relies on the `memcpy + SetSize` being compiled as-written and not as a shallow-copy-plus-compiler-optimization.  We believe the current IR cannot fold to a shallow copy, but this is not tested.  Mitigation: R-302 Part B (promote to MFEM_VERIFY).
- MFEM's partitioning of `BooleanFragments`-generated fault meshes — whether future MFEM versions or alternate mesh tools (e.g., reading a `.pmsh` with explicit shared-face tagging) will classify TPV102 faults as shared. If they do, R-001 activates silently, and if R-101/R-102 have bugs (R-301/R-305 documented above), production runs will silently drift.
- The four dev-queue sbatch scripts' actual Frontera environment — not verified on the machine; only the command-line invocation is reviewed.
