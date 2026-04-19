# Code Review v2: Fresh adversarial review of v1 fix round

Reviewing the code *after* `tpv102_debug_v1_fix.md` was applied.  This is
a full from-scratch audit — every pass re-executed, not a checklist of
prior findings.

## Review Scope
- Input docs:
  - `tpv102_debug_v1.md` (original bug report)
  - `tpv102_debug_v1_check.md` (first-round review, 7 findings)
  - `tpv102_debug_v1_fix.md` (first-round fix report)
- Files reviewed (diffed against `git HEAD`):
  - `miniapps/seas/drivers/tpv102_driver.cpp`
  - `miniapps/seas/dynamic/tpv102_setup.hpp`
  - `miniapps/seas/dynamic/wave_operator.hpp`
  - `miniapps/seas/dynamic/wave_operator.inl`
  - `miniapps/seas/io/paraview_output.hpp`
  - `miniapps/seas/.gitignore`
- Domain context consulted:
  - `miniapps/seas/CLAUDE.md` (sign-convention checklist)
  - `miniapps/seas/dynamic/fault_face_flux.cpp` (Evaluate symmetry properties)
  - MFEM public headers (`pmesh.hpp:460-580`) to verify
    `face_nbr_elements_offset` / `GetFaceNbrRank` are public
  - SCEC TPV102 spec (whole-space friction convention)

## Findings

### [R-101] [CRITICAL] [POSSIBLE] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — R-001's (+,−) swap assumes MFEM gives **identical** face normals on both ranks; the fix report flags this as unverified, and under the opposite convention the fix does not actually equalise `DOFData`

**Category:** ASSUMPTION

**Description:**
The fix swaps `(Q_plus_local, Q_minus_local)` on the non-owner rank so that `Evaluate` sees the same `(+,−)` arguments on both ranks. The claim in the fix report is that this makes rank A's and rank B's `DOFData` updates bit-identical. That claim is only true **if** `Tinv_A == Tinv_B` on both ranks — i.e., `CalcOrtho(ftr->Face->Jacobian())` produces the identical geometric normal regardless of which rank owns the face.

If MFEM re-orients the face so the normal always points "out of local Elem1" (the convention regular-interior faces rely on for `flux_.Interior` conservation), then `nor_A = −nor_B`, `Tinv_A ≠ Tinv_B`, and even with the swap:

```
Rank A: Evaluate(fdata_A, Tinv_A · Q_A, Tinv_A · Q_B, …)
Rank B: Evaluate(fdata_B, Tinv_B · Q_A, Tinv_B · Q_B, …)   (after swap)
```

The tangential-velocity components (`VY`, `VZ`) and the normal-velocity component (`VX`) pick up sign flips between `Tinv_A·Q` and `Tinv_B·Q`, so `ComputeTrialTraction`'s velocity-jump term (Eq. 7b) still differs, and the `DOFData` entries still diverge. Explicitly, on a planar fault with `nor_A = −nor_B`:

```
tau1_trial_A  ∝ (v_t1⁻  − v_t1⁺)   [in A's fault-local frame]
tau1_trial_B  ∝ (v_t1⁻' − v_t1⁺')  [in B's fault-local frame]

v_t1' = v · t1_B.  If BuildFrame gives t1_B = t1_A, the velocity-jump
term is the same; if BuildFrame gives t1_B = −t1_A (which is the
natural right-handed choice when nor flips and t2 is held fixed), the
term flips sign.
```

The fix report (§ "Unresolved Findings") explicitly lists this as an open question:

> "Orientation of `CalcOrtho(ftr->Face->Jacobian())` on shared faces.
>  R-001's fix is robust to either convention…"

Actually, **it is not**: the fix is correct only under the convention where `Tinv_A = Tinv_B`. The 4-rank 1000 m local test (`tpv102_debug_v1_fix.md § "Local Test Evidence"`) shows cross-rank `||Q||_∞` growth, which is *evidence* that the fix improved something, but does **not** prove `DOFData` on the two halves of each shared fault QP is consistent — it only shows that the two halves did not zero out wave propagation entirely. Drift of a few percent per step is consistent with what that test measured (r3 at ~1e-14, r1 at ~5e-12).

**The 400-rank 200 m production run is the actual failure mode we are trying to fix.** Running it and seeing whether the rupture propagates is the only dispositive check; the fix's correctness is a prerequisite for that run not wasting node-hours.

**Trigger:**
Any run of TPV102 where the hypocenter is separated from a far-field probe by at least one partition seam — i.e., every realistic production run.

**Actual behavior:**
`DOFData` for a shared fault QP *may* drift between the two ranks that own it, depending on MFEM's internal normal convention.  Fix report's claim of bit-identical updates is unverified.

**Expected behavior:**
`DOFData.{V1, V2, tau1_corr, tau2_corr, sigma_n_corr, psi}` at every shared fault QP must be bit-identical on both ranks that own it, on every RK4 stage, for every step, regardless of MFEM version.

**Suggested fix (two parts — must do both):**

**Part A — verify the MFEM invariant directly.** Add a one-shot diagnostic at the top of `ComputeSharedFaceFluxRHS` that dumps, for the first shared fault face on each rank, `nor` and `Tinv` + a signature of (local vertex indices, local coords), then `MPI_Allgather` it at step 0. Rank 0 asserts consistency. If the assertion fails, abort with a clear error telling the user "MFEM shared-face normal convention not what R-001 assumes — fix needed in `ComputeSharedFaceFluxRHS`." This is cheap (one step) and catches a future MFEM change immediately.

**Part B — add the DOFData consistency test as a dedicated MPI unit test.** See test_case below. Concretely:

```diff
 // In ComputeSharedFaceFluxRHS, after the sf loop, add:
+#ifdef MFEM_DEBUG
+   // One-shot consistency diagnostic: assert R-101 invariant at step 0.
+   static bool checked = false;
+   if (!checked && fault_dof_data_ && !shared_fault_dof_offset_.empty()) {
+      checked = true;
+      // Pick the first shared fault face on this rank and send its
+      // {psi, V1, tau1_corr} to the peer.  Peer does the same.  Both
+      // assert bit-equality.  (Implementation is ~30 lines; omitted for
+      // brevity here — see the test case below.)
+      verify_shared_fault_dofdata_consistency(*fault_dof_data_, pmesh);
+   }
+#endif
```

**Test case:**
```python
def test_R101_shared_fault_dofdata_bit_identical():
    """
    Minimal reproducer: 2-rank, 2-tet mesh with one shared fault face.
    Apply a nonzero Q perturbation on one side, run 1 RK4 stage, gather
    DOFData from both ranks via MPI, assert bit-equality.

    This directly verifies R-001 without the full TPV102 setup.
    """
    # mesh: two tets sharing a triangular face at y=0; partition so each
    # tet lives on its own rank.  Tag the shared face as fault (attr=3).
    mesh_path = build_two_tet_fault_mesh(tmp_path)

    # Initial state: nonzero Q on rank 0's element, zero on rank 1's.
    Q_init = perturbation_on_rank_zero_only()
    dof_data_ranks = run_one_rk4_stage(mesh_path, Q_init, nranks=2)

    for field in ["V1", "V2", "tau1_corr", "tau2_corr",
                  "sigma_n_corr", "psi"]:
        vals = gather_shared_fault_dofdata(dof_data_ranks, field)
        # Both ranks own the QP; values must be exactly equal.
        np.testing.assert_array_equal(
            vals[0], vals[1],
            err_msg=f"R-001 fix did NOT equalise DOFData.{field} — "
                    f"rank 0 has {vals[0]}, rank 1 has {vals[1]}")
```

---

### [R-102] [MODERATE] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — `peer_rank == −1` fallback silently masks resolution failure and re-introduces the R-001 bug

**Category:** BUG

**Description:**
```cpp
int peer_rank = (sf < static_cast<int>(shared_face_peer_rank_.size()))
              ? shared_face_peer_rank_[sf] : -1;
const bool owner = (peer_rank < 0) || (my_rank < peer_rank);
if (owner) { … /* no swap */ }
else       { … /* swap */ }
```

If peer-rank resolution fails in the constructor (`fn` loop doesn't find the face neighbor for whatever reason — MFEM API change, malformed mesh, etc.), `shared_face_peer_rank_[sf] = −1`. Here that `−1` is treated as "I'm the owner, no swap needed." If **both** ranks sharing the face hit the same resolution failure (likely, since the bug is symmetric), **both** think they are the owner, both skip the swap, and we are back to the exact pre-R-001 H3 situation — silently. No error, no warning.

The constructor already asserted the face neighbor was resolvable (by virtue of `ExchangeFaceNbrData` succeeding). A `−1` in `shared_face_peer_rank_` is a programmer-level bug that should abort the run, not fall back to buggy behaviour.

**Trigger:**
Any future change (MFEM upgrade, unusual mesh, weird partitioning) that causes the `fn` lookup loop in the ctor to fail for at least one shared face.

**Actual behavior:**
Silently falls back to unswapped Evaluate → shared-fault DOFData drift.

**Expected behavior:**
Hard abort with a descriptive error.

**Suggested fix:**
```diff
 int peer_rank = (sf < static_cast<int>(shared_face_peer_rank_.size()))
               ? shared_face_peer_rank_[sf] : -1;
-const bool owner = (peer_rank < 0) || (my_rank < peer_rank);
+MFEM_VERIFY(peer_rank >= 0,
+            "shared_face_peer_rank_[" << sf << "] was not resolved in "
+            "the WaveOperator ctor.  R-001 fix cannot canonicalise "
+            "the (+,−) side — aborting to avoid silently returning "
+            "the pre-fix buggy behaviour.");
+const bool owner = (my_rank < peer_rank);
```

Note: tie-break at `my_rank == peer_rank` cannot happen (a rank is never its own peer on a shared face), so the strict `<` is correct. For safety, also reject equality:

```diff
+MFEM_VERIFY(my_rank != peer_rank,
+            "shared face peer_rank equals my_rank — ParMesh invariant "
+            "violated");
```

**Test case:**
```python
def test_R102_unresolved_peer_rank_aborts():
    """
    Force a resolution failure by mocking GetFaceNbrRank to return -1.
    Expect MFEM_VERIFY abort, not silent buggy-path fallback.
    """
    # Since this is a sanity guard, test at unit level with a stub.
    with pytest.raises(RuntimeError, match="was not resolved"):
        run_wave_operator_with_stub_peer_rank_minus_one()
```

---

### [R-103] [MODERATE] `drivers/tpv102_driver.cpp` — Missing `#include <limits>` and `<vector>`; works only via transitive MFEM headers

**Category:** QUALITY (near-BUG)

**Description:**
The driver now uses `std::numeric_limits<real_t>::max()` (line 400) and `std::vector<int> watch`, `std::vector<int> watch_unique` (lines ~820-830) but does not `#include <limits>` or `<vector>` directly. It compiles today only because `mfem.hpp` transitively includes both. Any MFEM header cleanup that drops transitive `<limits>` or `<vector>` breaks this driver silently (missing-symbol link error at best, worse compile error).

This is not esoteric: CI systems with aggressive header-include-what-you-use linters (e.g., `iwyu`) flag exactly this.

**Trigger:**
Any future MFEM header refactor that tightens includes; building with `-ftrack-macro-expansion=0` or strict IWYU.

**Actual behavior:**
Compiles only by luck.

**Expected behavior:**
Direct `#include` of every stdlib header the TU uses.

**Suggested fix:**
```diff
 #include <algorithm>
 #include <cstring>
 #include <iostream>
 #include <fstream>
+#include <limits>
 #include <memory>
 #include <string>
 #include <cmath>
+#include <vector>
 #include <sys/stat.h>
```

**Test case:** N/A — verified by compile step under `-fsanitize=leak` with headers forcibly reduced.  Track via a one-line lint check in CI:

```yaml
- name: IWYU-lite
  run: |
    grep -n 'std::numeric_limits\|std::vector' miniapps/seas/drivers/*.cpp |
      while read loc; do
        f=${loc%%:*}
        grep -q '^#include <limits>' $f || { echo "$f: missing <limits>"; exit 1; }
      done
```

---

### [R-104] [MODERATE] `drivers/tpv102_driver.cpp:paraview_write` — `PeekShouldWrite` + `Save`/`ShouldWrite` evaluate the schedule twice with *possibly inconsistent* arguments; drift is possible under `OutputInterval(V_max)` path

**Category:** BUG (edge case)

**Description:**
The paraview_write lambda now evaluates the schedule up to three times per step (Peek, then Save *or* ShouldWrite). For step-based scheduling this is harmless — `cycle % N == 0` is deterministic. For **time-based** scheduling via `OutputInterval(V_max)` (the adaptive V-keyed schedule used when neither `--paraview-every` nor `--paraview-dt` is set), `OutputInterval` picks between 0.01 s (coseismic), 1 s (nucleation), and 1 year (interseismic) buckets based on `V_max` — a **step-function** of V_max.

If `V_max` straddles a bucket boundary across the two calls (possible if anything mutates shared state between Peek and Save/ShouldWrite — e.g., the fault packing loop — but in this code nothing does), Peek might return true using the coseismic 0.01 s interval while Save uses the nucleation 1 s interval, or vice versa, deciding inconsistently whether `time - last_write_time_` clears the gate. Even if they agree on the same V_max (which they should, by inspection), this introduces a new invariant that the next maintainer must preserve — otherwise the ParaView output drops cycles silently.

**Trigger:**
Time-based scheduling (default when `use_paraview` is on and neither interval flag is set). If anyone ever passes *different* V_max values to Peek and Save — or if `last_write_time_` gets mutated between them — the schedule diverges.

**Actual behavior:**
Schedule gate evaluated twice. Correct today by accident (no shared-state mutation between the two calls, identical V_max), but fragile.

**Expected behavior:**
Evaluate the schedule **once** per `paraview_write` call.

**Suggested fix — collapse Peek and Save into one call that both gates and records:**

```diff
 auto paraview_write = [&](int step_num, real_t time, real_t V_max)
 {
    if (!pv_out) { return; }

-   // R-007 fix: skip packing on cycles that are not scheduled frames.
-   // PeekShouldWrite is idempotent — does NOT advance last_write_time_.
-   // The subsequent Save()/ShouldWrite() call still performs the write
-   // and updates the schedule.
-   if (!pv_out->PeekShouldWrite(step_num, time, V_max)) { return; }
-
-   // When volume PVD is suppressed, we still need to decide whether
-   // this cycle is a scheduled frame (for fault-surface VTUs).
-   bool wrote = false;
-   if (pv_no_domain)
-   {
-      wrote = pv_out->ShouldWrite(step_num, time, V_max);
-   }
-   else
-   {
-      std::memcpy(pv_vel_gf->GetData(), …);
-   }
-   // … packing …
-   if (!pv_no_domain) {
-      pv_out->UpdateFaultFieldsBP5(…);
-      wrote = pv_out->Save(step_num, time, V_max);
-   }
+   // Single-shot schedule gate:
+   if (!pv_out->PeekShouldWrite(step_num, time, V_max)) { return; }
+
+   // Always pack fault fields (cheap; needed for fault-surface VTU).
+   // Copy velocity only when writing the volume PVD.
+   if (!pv_no_domain) {
+      std::memcpy(pv_vel_gf->GetData(), …);
+   }
+   pack_fault_fields();   // extract the per-QP loop into a helper
+
+   // Commit the write and advance the schedule.
+   bool wrote = pv_no_domain
+              ? pv_out->CommitSchedule(step_num, time, V_max)  // advances last_write_time_, returns true
+              : pv_out->Save(step_num, time, V_max);
+   // (CommitSchedule is a new helper: same body as ShouldWrite but
+   // always returns true.  Trivially defined next to PeekShouldWrite.)
```

Add to `paraview_output.hpp`:
```diff
+   /// Advance last_write_time_ to `time` and return true.
+   /// Use after PeekShouldWrite has already confirmed the cycle is
+   /// scheduled.  Idempotent under repeat calls with the same (cycle,time).
+   bool CommitSchedule(int cycle, real_t time, real_t /*V_max*/)
+   {
+      last_write_time_ = time;
+      return true;
+   }
```

**Test case:**
```python
def test_R104_paraview_schedule_evaluated_once():
    """
    Run with adaptive V_max scheduling and verify the schedule is
    evaluated exactly once per paraview_write — no double-eval.
    Instrument ParaViewOutput with a call counter on PeekShouldWrite
    and Save/ShouldWrite; assert sum equals nsteps.
    """
    counter = InstrumentedPVOutput()
    run_driver_with_counter(counter, nsteps=100, adaptive_schedule=True)
    # Expect exactly 1 schedule evaluation per step:
    assert counter.peek_calls + counter.save_calls == counter.total_steps
```

---

### [R-105] [MODERATE] `drivers/tpv102_driver.cpp:main` — `hypo_rank` detection relies on an MPI struct layout that is not guaranteed to match `MPI_DOUBLE_INT`

**Category:** ASSUMPTION

**Description:**
```cpp
struct { double d; int r; } in{local_min_dist2, rank}, out{};
MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
```

`MPI_DOUBLE_INT` is defined by the MPI standard as the layout of a specific struct defined in `mpi.h`, typically `struct { double d; int i; }`. C++ does not guarantee that `struct { double d; int r; }` (your anonymous struct) has the same alignment and padding as `mpi.h`'s version — in practice it does on every compiler we use, but it's a classic portability footgun. The `struct{...}` will almost certainly be laid out as `{double, int, 4-byte-pad}` = 16 bytes, same as MPI_DOUBLE_INT, but nothing in the language standard forces that.

More important: `real_t` may be `float` in a mixed-precision build (`MFEM_USE_SINGLE`). The line `in{local_min_dist2, rank}` silently widens a `float` to `double` — that's fine as an implicit conversion — but the driver elsewhere uses `real_t` consistently, and this spot quietly breaks that pattern. If someone refactors to `struct { real_t d; int r; }` + `MPI_DOUBLE_INT`, it silently mis-decodes on single-precision builds.

**Trigger:**
Future build with `MFEM_USE_SINGLE=YES`; or compiler/platform with unusual `{double,int}` struct padding.

**Actual behavior:**
Depends on compiler layout of an anonymous struct against the MPI_DOUBLE_INT layout hardcoded in `mpi.h`.

**Expected behavior:**
Use `std::pair<double,int>` or a named struct with explicit `static_assert(sizeof(S) == sizeof(double) + sizeof(int) + padding)`; or use MPI's provided type `struct { double; int; }` from the MPI header directly (not all MPI implementations expose this).

**Suggested fix:**
```diff
+   struct MinDist { double d; int r; };
+   static_assert(offsetof(MinDist, d) == 0, "MinDist layout != MPI_DOUBLE_INT");
+   static_assert(offsetof(MinDist, r) == sizeof(double),
+                 "MinDist layout != MPI_DOUBLE_INT (unexpected padding)");
-   struct { double d; int r; } in{local_min_dist2, rank}, out{};
+   MinDist in{static_cast<double>(local_min_dist2), rank};
+   MinDist out{};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    hypo_rank = out.r;
```

**Test case:** compile-time only; the `static_assert`s catch the bug at build time.

---

### [R-106] [LOW] `drivers/tpv102_driver.cpp:debug-qnorm watch list` — When `hypo_rank` is near `nprocs−1`, the watch list collapses and the diagnostic shows only one rank

**Category:** QUALITY

**Description:**
```cpp
std::vector<int> watch = {
   hypo_rank,
   std::min(hypo_rank + 1, nprocs - 1),
   std::min(hypo_rank + 4, nprocs - 1),
   nprocs - 1};
```

For `hypo_rank == nprocs − 1`, all four entries collapse to `nprocs − 1` and the dedup'd watch list has a single rank — the output line becomes just `[qnorm:watch] r399=… (hypo_rank=399)` and the user has no neighbours to compare against. The diagnostic's stated purpose ("H1 diagnosis") is lost precisely when `hypo_rank` happens to be at the high end.

For TPV102 specifically this is unlikely (hypocenter is near the geometric centre, METIS usually places it mid-range), but it's not impossible.

**Trigger:**
Partition that places the hypocenter at the last rank.

**Actual behavior:**
Diagnostic collapses to a single-rank watch.

**Expected behavior:**
Always show at least 2–3 ranks from the hypocenter neighborhood plus a far rank.

**Suggested fix:**
```diff
-std::vector<int> watch = {
-   hypo_rank,
-   std::min(hypo_rank + 1, nprocs - 1),
-   std::min(hypo_rank + 4, nprocs - 1),
-   nprocs - 1};
+std::vector<int> watch = { hypo_rank };
+for (int off : {1, 4})                     // neighbours "above"
+{
+   if (hypo_rank + off <  nprocs) watch.push_back(hypo_rank + off);
+   if (hypo_rank - off >= 0)      watch.push_back(hypo_rank - off);
+}
+watch.push_back(nprocs - 1);               // far rank
+watch.push_back(0);                        // always include rank 0
```

**Test case:** Unit-test the watch-list construction against edge cases:
```python
@pytest.mark.parametrize("nprocs,hypo,expected", [
    (1,   0,   [0]),
    (4,   0,   [0, 1, 2, 3]),
    (4,   3,   [3, 2, 1, 0]),             # all ranks visible even at top
    (400, 399, [399, 398, 395, 0]),       # neighbours to the left, plus r0
])
def test_R106_watch_list_non_collapsing(nprocs, hypo, expected):
    assert sorted(build_watch_list(nprocs, hypo)) == sorted(expected)
```

---

### [R-107] [LOW] `dynamic/wave_operator.inl:ctor` — `shared_face_peer_rank_.assign(n_shared, -1)` happens **before** `GetSharedFaceTransformations(sf)` returns nullptr; the −1 sentinel coexists with two meanings (unresolved / not-a-fault) and R-102's check cannot distinguish them

**Category:** QUALITY

**Description:**
`shared_face_peer_rank_[sf]` is initialized to `-1` for every shared face.  It is **only overwritten** when `ftr = GetSharedFaceTransformations(sf)` is non-null **and** `fn` lookup succeeds.  For shared faces whose `ftr` is null (can happen for "ghost-only" shared faces in certain MFEM configurations), the value stays at `-1` even though the face is not a fault face and the fix never consults it. That is fine today.

But R-102's suggested abort (`MFEM_VERIFY(peer_rank >= 0, …)`) only fires for *fault* shared faces — so R-102 works. What breaks R-107 is future code that consults `shared_face_peer_rank_[sf]` for a non-fault face and misinterprets the sentinel. Separate the two conditions now with a distinct sentinel or a `std::optional<int>`:

**Suggested fix:**
Either (a) use `std::optional<int>`, or (b) split into `bool ftr_valid_` + `int peer_rank_`:

```diff
-std::vector<int> shared_face_peer_rank_;
+struct SharedFacePeer {
+   bool  ftr_valid = false;   // true iff GetSharedFaceTransformations succeeded
+   int   peer_rank = -1;      // valid only when ftr_valid
+};
+std::vector<SharedFacePeer> shared_face_peer_;
```

**Test case:** covered by R-102's test.

---

### [R-108] [LOW] `io/paraview_output.hpp:PeekShouldWrite` — signature is `const` but `Save`/`ShouldWrite` are **not** const; subtle API mismatch

**Category:** QUALITY

**Description:**
`PeekShouldWrite` is declared `const` (correctly — it doesn't mutate).  `Save` and `ShouldWrite` aren't.  Users reasonably assume the three are interchangeable gate-checks; the `const`/non-`const` split is a subtle trap for the next person who reads the code.  Either make `Peek…`/`ShouldWrite` both `const` and factor out the mutation into a separate `Advance` call (R-104's suggestion does this), or leave `Peek…` non-const to match (less desirable — wastes the `const` annotation).

**Suggested fix:** adopt R-104's `CommitSchedule` helper so the split is explicit (Peek = const read, Commit = non-const write).

**Test case:** compile-time only (const correctness).

---

### [R-109] [LOW] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — `my_rank` redeclared every call via `pmesh.GetMyRank()`; cache it in the operator

**Category:** QUALITY

**Description:**
Every call to `ComputeSharedFaceFluxRHS` (= every RK4 stage × every step — 4 × nsteps times) re-queries `pmesh.GetMyRank()`.  Cheap but wasteful; also, caching it in the ctor makes the value explicit as an invariant of the operator.  Drop into a member.

**Suggested fix:**
```diff
 // wave_operator.hpp:
+int my_rank_ = 0;

 // wave_operator.inl ctor:
+#ifdef MFEM_USE_MPI
+   if constexpr (IsParallelMesh<MeshType>::value)
+   {
+      my_rank_ = static_cast<const ParMesh &>(mesh_).GetMyRank();
+   }
+#endif

 // ComputeSharedFaceFluxRHS:
-int my_rank = pmesh.GetMyRank();
```

**Test case:** N/A — pure refactor.

---

## Recommended additional unit tests

These are the gaps **not filled** by the existing suite, ranked by value:

1. **R-101 direct regression test (HIGHEST PRIORITY).**  2-rank MPI test with a 2-tet fault mesh. Assert `DOFData` bit-equality on both ranks after 1 RK4 stage and after 10 RK4 steps.  **Not implementing this leaves R-001's fix unverified.**
2. **Shared-face normal convention assertion test.**  1-step MPI test that dumps `nor_vec` and `Tinv` for each shared fault face on both ranks; `MPI_Allgather` and cross-check. Use this both as a unit test and as a once-per-run startup diagnostic (see R-101 Part A suggestion).
3. **Godunov antisymmetry at a shared non-fault face.**  Existing `seas_test_parallel_wave_operator` uses a uniform cube; it does not exercise the *graded-mesh irregular partition* that TPV102 produces. Add a 4-rank test with a graded tet mesh (h varies by 10×) and an initial point-source pulse; assert that `||Q||_∞` on a rank far from the pulse grows above noise floor by t = h_min / cs × 10.
4. **Frontera coarse regression (400-rank, but 1000 m mesh).**  Even at `h = 1000 m` with `Λ_dyn / h ≪ 1` rupture will not propagate, but the *cross-rank bulk wave energy* diagnostic can be checked: use `--debug-qnorm` to verify ranks geographically distant from the hypocenter register nonzero `||Q||_∞` within tens of ms after nucleation engages. Cheap (~30 min wall) and decisive for H1.
5. **Shared-fault-face ParaView smoke.** Run 8-rank 1000 m TPV102 for 0.5 s with `--paraview --paraview-every 10`; open the fault-surface PVD in ParaView and check visually for **stripes of opposite sign** in `slip_rate_strike` at partition seams. If R-001 is correct, no stripes.  If R-001 is incomplete, stripes are the flagship visual signature.

## Recommended short Frontera dev runs

Dev queue = 2 h wall, no charge beyond node-hours.  In priority order:

1. **`devq-4rank-1000m-debug-qnorm` (15 min).**  4 ranks, `tpv102_1000m.msh`, `--order 1`, `tfinal=1.5 s`, `--debug-qnorm`.  Expected: `[qnorm:watch]` line shows non-hypo ranks at ≥ 1e-18 by t = 0.2 s (wave transits 1 km at `c_s = 3464 m/s`).  **Catches any residual R-101 bug without a full production run.**
2. **`devq-50rank-200m-debug-qnorm` (2 h, dev queue).**  50 ranks, `tpv102_200m.msh`, `--order 1`, `--dt` matching the CFL-auto result, `tfinal=1.5 s`, `--debug-qnorm --no-domain-pv --paraview --paraview-every 50`.  Expected: non-hypo ranks register `||Q||_∞ > 0` and rupture starts to extend past the nucleation patch.  Disk: fault-surface PVD only, ~50 MB.
3. **`devq-400rank-200m-repro-v1` (2 h, dev queue).**  Repeat the `results_200m_p1_3s_dev_job7664258` configuration with the v1 fix landed, with `--debug-qnorm`.  Expected: fault-surface PVD at t=1.5 s shows rupture at |x| ≥ 1.5 km — i.e., the v1 symptom is gone.  This is the dispositive test.  **Do not burn the production-queue node-hours until 1+2 pass.**
4. **Follow-up production run once 1–3 pass.**  200 m, 400 rank, `--order 1`, 3 s, production queue.  If the dev-queue runs all show cross-rank propagation, proceed to the 12 s production run.

## Summary
- Critical issues: 1 (R-101 — the R-001 fix is only partly verified; has a load-bearing unverified MFEM invariant)
- Moderate issues: 4 (R-102 silent peer-rank fallback, R-103 missing includes, R-104 double schedule evaluation, R-105 MPI_DOUBLE_INT layout assumption)
- Low issues: 4 (R-106 watch-list collapse, R-107 sentinel overload, R-108 const mismatch, R-109 my_rank caching)
- Plan compliance: **PARTIAL** — the 7 v1 findings are addressed, but R-001's fix is not directly tested, and R-005's fix is guarded only by a debug-mode `MFEM_ASSERT` that does not fire in Release builds.
- Verdict: **PASS WITH FIXES.** Must fix R-101 (add a direct test) and R-102 (hard-abort fallback) before the next Frontera production run.  R-103/R-104/R-105 are small and can ride along.

## Unreviewed Areas
- `--pv-low-order` interaction with `SetHighOrderOutput(false) + SetLevelsOfDetail(1)` in MFEM — we confirmed visually that the output size drops, but we have not reviewed whether MFEM's linear-output path projects the order-p velocity field correctly (vs. just dropping it to the P1 corner DOFs of each tet). This affects the **faithfulness** of ParaView velocity, not correctness of the simulation.  Out of scope for this review.
- Regression against Tandem / SeisSol TPV102 reference traces (SCEC cross-code).  Only applicable once the 400-rank production run completes without the H1 symptom.
- The new `tpv102/mesh/tpv102_1000m.geo` file (untracked).  Not reviewed — it's a new mesh input, not code.
