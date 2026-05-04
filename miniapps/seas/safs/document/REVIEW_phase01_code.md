# Code Review: Phase 0 + Phase 1 implementation (2026-04-29)

## Review Scope

- **Plan reviewed against:** `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/safs/PLAN_cgal_corefine.md` (Phase 0 + Phase 1)
- **Files reviewed:**
  - `tools/CMakeLists.txt` (143 lines)
  - `tools/main.cpp` (53 lines — Phase 0 hello binary)
  - `tools/README.md` (167 lines)
  - `tools/cgal_patch/README.md` (73 lines)
  - `tools/corefine_faults/corefine.hpp` (67 lines)
  - `tools/corefine_faults/corefine.cpp` (177 lines)
  - `tools/corefine_faults/io.hpp` (88 lines)
  - `tools/corefine_faults/io.cpp` (291 lines)
  - `tools/corefine_faults/main.cpp` (422 lines)
  - `tools/tests/test_corefine_smoke.cpp` (215 lines)
- **Domain context:** `seas-mfem/CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  `REVIEW_plan_cgal_corefine.md` (the plan-review the implementation
  was meant to satisfy), CGAL 5.6.1 PMP docs (corefine,
  isotropic_remeshing, repair_polygon_soup).
- **Live verification performed:**
  - Both targets build cleanly (`corefine_faults_hello`,
    `corefine_faults`, `test_corefine_smoke`).
  - `corefine_faults_hello` exits 0, prints `CGAL version: 5.6.1`.
  - `test_corefine_smoke` runs all 4 sub-tests and exits 0.
  - End-to-end run on Mill Creek × SBMT-SAF 2000 m fixture
    (`output/two_crossing_2000m/stl_raw/`):
    - 2 faults read (793 + 1346 triangles), 1 pair, runtime 21 ms.
    - `n_constrained_edges_A == n_constrained_edges_B == 143`.
    - Post-corefine triangle counts: 1079 + 1633.
    - `triangle_to_fault.json` and `intersection_report.json` emit
      with the documented schema.

I assumed at least 3 bugs and found 11 actionable findings (2
critical, 5 moderate, 4 low). Build infrastructure and the C++
recipe itself are sound; the issues are concentrated in (a) the
plan-vs-actual constrained-edge count mismatch, (b) a tautological
schema check, (c) silent failure paths in the gate logic, and (d)
a few small CGAL-API hygiene issues.

---

## Findings

### [R-001] [CRITICAL] [tools/corefine_faults/main.cpp:402–413] — `schema_validation` block is a tautology that always passes

**Category:** BUG

**Description:**
After writing the JSON sidecars, `main.cpp` runs what claims to
be a schema-invariant check:

```cpp
// 7. Schema-invariant check on triangle_to_fault.json.
{
    std::int64_t total = 0;
    for (const auto& f : per_fault) total += f.n_triangles;
    std::int64_t expected = 0;
    for (const auto& f : per_fault) expected += f.n_triangles;
    if (total != expected) {
        std::cerr << "[corefine_faults] triangle_to_fault schema "
                     "invariant violated.\n";
        return 3;
    }
}
```

`total` and `expected` are computed from the **same data via the
same loop**. They are bitwise-identical by construction. The
`if (total != expected)` branch is unreachable. Exit code 3
(schema-validation error) cannot fire.

The intent (per plan §Phase 1 §6 + §Acceptance criteria) is to
verify that the on-disk JSON's `n_total_triangles` equals the
sum of `faults[*].n_triangles` AND that ranges are contiguous.
The current code does neither — it just sums the in-memory
`per_fault` vector against itself.

**Trigger:** Always. The check is structurally a no-op.

**Actual behavior:** Exit code 3 is allocated in the plan but
never returned, even when the JSON writer corrupts the schema
(e.g., a future bug where `write_triangle_to_fault_json` skips
a fault entry).

**Expected behavior:** Re-read the just-written JSON, parse out
`n_total_triangles` and the `faults[*].n_triangles` values,
verify the sum equality AND verify the per-fault ranges form a
contiguous partition of `[0, n_total_triangles)`.

**Suggested fix:**

Add a small JSON re-reader (or use a hand-rolled scan that pulls
the integer fields by string match). Concretely, replace the
tautology with two real checks: one against the on-disk JSON,
one against the in-memory vector.

```diff
@@ tools/corefine_faults/main.cpp lines 402–413
-    // 7. Schema-invariant check on triangle_to_fault.json.
-    {
-        std::int64_t total = 0;
-        for (const auto& f : per_fault) total += f.n_triangles;
-        std::int64_t expected = 0;
-        for (const auto& f : per_fault) expected += f.n_triangles;
-        if (total != expected) {
-            std::cerr << "[corefine_faults] triangle_to_fault schema "
-                         "invariant violated.\n";
-            return 3;
-        }
-    }
+    // 7. Schema-invariant check on triangle_to_fault.json (R-001):
+    //    re-read the just-written JSON and verify (a) sum of
+    //    per-fault n_triangles equals n_total_triangles, and
+    //    (b) ranges form a contiguous partition.
+    {
+        std::ifstream fh(tri2fault);
+        if (!fh) {
+            std::cerr << "[corefine_faults] cannot re-open "
+                      << tri2fault << " for schema check.\n";
+            return 3;
+        }
+        const std::string body((std::istreambuf_iterator<char>(fh)),
+                                std::istreambuf_iterator<char>());
+        // Pull the n_total_triangles value with a strict scan.
+        const auto pos = body.find("\"n_total_triangles\":");
+        if (pos == std::string::npos) {
+            std::cerr << "[corefine_faults] schema: "
+                         "n_total_triangles key missing.\n";
+            return 3;
+        }
+        std::int64_t on_disk_total = 0;
+        if (std::sscanf(body.c_str() + pos,
+                         "\"n_total_triangles\": %lld",
+                         reinterpret_cast<long long*>(&on_disk_total)) != 1) {
+            std::cerr << "[corefine_faults] schema: "
+                         "n_total_triangles parse failed.\n";
+            return 3;
+        }
+        std::int64_t in_mem_sum = 0;
+        for (const auto& f : per_fault) in_mem_sum += f.n_triangles;
+        if (on_disk_total != in_mem_sum) {
+            std::cerr << "[corefine_faults] schema: on-disk total "
+                      << on_disk_total << " != in-memory sum "
+                      << in_mem_sum << '\n';
+            return 3;
+        }
+        // Range contiguity: each fault's range[0] must equal the
+        // previous fault's range[1], starting at 0.
+        std::int64_t cursor = 0;
+        for (const auto& f : per_fault) {
+            // Find this fault's range entry; verify range == [cursor, cursor + n_triangles].
+            const auto fkey = "\"" + f.short_name + "\":";
+            const auto fpos = body.find(fkey);
+            if (fpos == std::string::npos) {
+                std::cerr << "[corefine_faults] schema: fault "
+                          << f.short_name << " missing from JSON.\n";
+                return 3;
+            }
+            std::int64_t lo = -1, hi = -1;
+            if (std::sscanf(body.c_str() + fpos,
+                "\"%*[^\"]\": {\"n_triangles\": %*lld, \"range\": [%lld, %lld]}",
+                reinterpret_cast<long long*>(&lo),
+                reinterpret_cast<long long*>(&hi)) != 2) {
+                std::cerr << "[corefine_faults] schema: fault "
+                          << f.short_name << " range parse failed.\n";
+                return 3;
+            }
+            if (lo != cursor || hi != cursor + f.n_triangles) {
+                std::cerr << "[corefine_faults] schema: fault "
+                          << f.short_name << " range [" << lo << ", "
+                          << hi << "] does not match expected ["
+                          << cursor << ", " << cursor + f.n_triangles
+                          << "].\n";
+                return 3;
+            }
+            cursor = hi;
+        }
+        if (cursor != on_disk_total) {
+            std::cerr << "[corefine_faults] schema: cumulative "
+                      << cursor << " != on-disk total "
+                      << on_disk_total << '\n';
+            return 3;
+        }
+    }
```

If the implementer prefers not to reparse JSON, accept that the
on-disk check is slightly more involved than the in-memory loop
shows; the alternative is to bake range computation into a
single function and verify by computing it twice in two
DIFFERENT ways. The current "two identical loops" formulation
is structurally inert.

**Test case:**
```cpp
// tools/tests/test_schema_invariant.cpp (new file)
TEST(R001) {
  // Run main with two faults, then mutate the on-disk JSON to
  // drop one fault's range entry.  Re-run an embedded validator
  // and assert it catches the corruption.  Confirms the check
  // can actually fire.
  // ... (smoke setup omitted for brevity) ...
  const fs::path j = out_dir / "triangle_to_fault.json";
  std::string body = read_file(j);
  // Corrupt: replace n_total_triangles with a wrong value.
  body = std::regex_replace(body,
      std::regex("\"n_total_triangles\": \\d+"),
      "\"n_total_triangles\": 999999");
  write_file(j, body);
  // Run only the validator (factor it out into a function).
  ASSERT_EQ(validate_triangle_to_fault_json(j, per_fault),
            /*exit code*/ 3);
}
```

---

### [R-002] [CRITICAL] [PLAN_cgal_corefine.md vs actual run] — Plan acceptance "n_constrained_edges ≥ 500 floor" is not met by the implementation (actual: 143)

**Category:** DEVIATION (plan-vs-implementation mismatch — needs
either plan update or under-count investigation)

**Description:**
PLAN_cgal_corefine.md Phase 1 acceptance criteria states:

> On the Mill Creek × SBMT-SAF 2000 m fixture
> (`output/two_crossing_2000m/stl_raw/`), the tool runs in ≤ 30 s
> on a laptop and reports `n_constrained_edges_A ==
> n_constrained_edges_B`, both > 200 (ballpark: the Python
> pipeline reports a single polyline of 581 vertices → ~580
> edges, so post-corefine constrained-edge count ≥ 500 is the
> floor).

The actual end-to-end run (verified against
`output/two_crossing_2000m/stl_raw/`):

```
[corefine_faults] pair safs_sbmt_millcreek x safs_sbmt_saf:
        ce_A=143 ce_B=143 (pre A,B=793,1346 → post=1079,1633)
```

143 is below the plan's stated `> 200` and well below the
`≥ 500 floor`. The criterion is unmet.

**Why the plan estimate was off (most likely):**
- Python's `chain_segments` produces polyline VERTICES, not the
  geometric polyline EDGES; with snap-key bucketing at 1 cm the
  Python pipeline recorded 581 polyline vertices, but many of
  those are sub-edge subdivisions that CGAL's exact corefine
  collapses into single edges.
- The plan extrapolated 580 polyline-edges (= 581 vertices − 1)
  from the Python output, treating polyline-vertex count as a
  one-to-one proxy for constrained-edge count. CGAL's measure
  is the count of intersection EDGES in the corefined
  triangulation; these are by-construction the topologically
  necessary edges and should be FEWER than Python's snap-induced
  vertices.

So 143 is plausibly **the correct geometric answer** and the
plan's 500-floor was based on a mismatched proxy. But this needs
to be VERIFIED by direct measurement (e.g., compute polyline
total length from CGAL output vs from Python output; if they
agree, CGAL is correct and Python over-segments).

**Trigger:** The current real-fixture run.

**Actual behavior:** The C++ smoke test passes (only checks
`> 0`), but the plan acceptance fails. There is no automated
gate that would catch this on CI today; only the plan document
disagrees with the implementation.

**Expected behavior:** Either (a) the plan is updated with the
empirically-measured count + a justification, OR (b) the
implementer demonstrates that 143 is wrong and CGAL is
under-counting.

**Suggested fix (option A — preferred unless evidence
suggests CGAL is buggy):**

```diff
@@ PLAN_cgal_corefine.md Phase 1 acceptance criteria
- On the Mill Creek × SBMT-SAF 2000 m fixture
   (`output/two_crossing_2000m/stl_raw/`), the tool runs in ≤ 30 s
   on a laptop and reports
   `n_constrained_edges_A == n_constrained_edges_B`, both > 200
-  (ballpark: the Python pipeline reports a single polyline of
-  581 vertices → ~580 edges, so post-corefine constrained-edge
-  count ≥ 500 is the floor).
+  (CGAL EPICK on this fixture produces 143 constrained edges
+  per side; the Python pipeline's "581 polyline vertices" is
+  an over-count caused by Python's snap-induced sub-edges.
+  The geometrically meaningful invariant is symmetry —
+  `n_constrained_edges_A == n_constrained_edges_B` — plus
+  a numeric floor of `> 50` to detect the case where corefine
+  finds nothing.  Empirical pin: 143 ± 10 on the 2000 m
+  fixture; record this in `intersection_report.json` so future
+  CGAL upgrades can be A/B compared.).
```

**Suggested fix (option B — if option A is rejected):**

Add a debug flag that emits the polyline as a DXF or VTK file,
compute its total Euclidean length, compare against an
independent measurement of the Python pipeline's polyline
length. If the lengths agree to within 1%, option A is correct.
If CGAL's polyline is materially shorter, there's an
under-count bug and corefine is missing intersection segments.

**Test case (option A):**
```cpp
TEST(R002) {
  // Pin the empirically-measured count so future CGAL upgrades
  // surface as a gate failure rather than silent drift.
  const auto cr = run_on_real_fixture(/*mill_creek*/, /*sbmt_saf*/);
  EXPECT_EQ(cr.n_constrained_edges_A, cr.n_constrained_edges_B);
  EXPECT_GT(cr.n_constrained_edges_A, 50);
  EXPECT_NEAR(cr.n_constrained_edges_A, 143, 10);
}
```

---

### [R-003] [MODERATE] [tools/corefine_faults/main.cpp:336–337 + manifold_gate semantics] — Manifold-gate FAIL string is recorded but never causes a non-zero exit

**Category:** BUG (silent failure)

**Description:**
`manifold_gate` (lines 178–186 of `main.cpp`) returns one of
three strings:

```cpp
std::string manifold_gate(const Mesh& m) {
    if (!CGAL::is_valid_polygon_mesh(m)) {
        return "FAIL_INVALID_POLYGON_MESH";
    }
    if (PMP::does_self_intersect(m)) {
        return "FAIL_SELF_INTERSECT";
    }
    return "PASS";
}
```

In the main loop (lines 336–337):
```cpp
ps.gate_manifold_A = manifold_gate(meshes[i]);
ps.gate_manifold_B = manifold_gate(meshes[j]);
```

Both gates are stored in the `PairStats` struct, written to
`intersection_report.json`, and otherwise **ignored**. There is
no `if (ps.gate_manifold_A != "PASS") return 2;` block. The
binary exits 0 even with FAIL strings in the JSON.

This is a regression vs the plan, which (per §Phase 1 §6 +
P-008 fixes from `REVIEW_plan_cgal_corefine.md`) treats
non-manifold output as a fatal corefine-runtime error (exit 2).

**Trigger:** Any input that drives `corefine` to produce a
post-mutation mesh where `is_valid_polygon_mesh` returns false
(possible if input has near-coincident triangles that corefine
mishandles, even if `throw_on_self_intersection` didn't fire on
the input). With current Mill Creek × SBMT-SAF the gate is PASS,
but the failure path is silently ignored.

**Actual behavior:** JSON records `"manifold_A": "FAIL_..."`;
exit code is 0; downstream `generate_safs_mesh.py` sees a
working `triangle_to_fault.json` and proceeds, possibly handing
HXT a non-manifold STL.

**Expected behavior:** Any FAIL string in either gate is fatal;
exit 2 with a message naming the offending fault and gate.

**Suggested fix:**

```diff
@@ tools/corefine_faults/main.cpp, after line 337
             ps.gate_manifold_A = manifold_gate(meshes[i]);
             ps.gate_manifold_B = manifold_gate(meshes[j]);
+            if (ps.gate_manifold_A != "PASS") {
+                std::cerr << "[corefine_faults] manifold gate FAILED "
+                          << "on fault " << ps.short_a
+                          << " after corefine with " << ps.short_b
+                          << ": " << ps.gate_manifold_A << '\n';
+                return 2;
+            }
+            if (ps.gate_manifold_B != "PASS") {
+                std::cerr << "[corefine_faults] manifold gate FAILED "
+                          << "on fault " << ps.short_b
+                          << " after corefine with " << ps.short_a
+                          << ": " << ps.gate_manifold_B << '\n';
+                return 2;
+            }
```

And similarly for `gate_interior_only` — currently set to
`"FAIL"` on idempotency violation but not propagated:

```diff
@@ tools/corefine_faults/main.cpp, after line 347
             ps.gate_interior_only = interior_crossing_only(
                 meshes[i], meshes[j], ps.short_a, ps.short_b)
                 ? "PASS" : "FAIL";
+            if (ps.gate_interior_only != "PASS") {
+                std::cerr << "[corefine_faults] interior-crossing-"
+                             "only gate FAILED on pair (" << ps.short_a
+                          << ", " << ps.short_b
+                          << "): re-corefine introduced new triangles, "
+                             "indicating the first corefine did not "
+                             "produce a fully-conformal pair.\n";
+                return 2;
+            }
```

**Test case:**
```cpp
TEST(R003) {
  // Inject a synthetic mesh where corefine's output happens to
  // have a self-intersection (e.g. a deliberately-crafted
  // near-coplanar pair).  Run main(); assert exit code 2.
  ASSERT_EQ(run_main(...), 2);
  // And verify the JSON still has the FAIL string for diagnostic.
  auto rep = parse_report(...);
  EXPECT_EQ(rep["pairs"]["A__x__B"]["gates"]["manifold_A"],
            "FAIL_SELF_INTERSECT");
}
```

---

### [R-004] [MODERATE] [tools/corefine_faults/main.cpp:339–342 + smoke test] — `polyline_edge_coincidence` gate hard-codes "PASS" by tautology; never actually verified

**Category:** ASSUMPTION (untested invariant)

**Description:**
Lines 339–342 of `main.cpp`:
```cpp
// Polyline-edge-coincidence gate is a tautology of corefine
// (the constrained edges are by construction the same on
// both sides of the cut), so we record PASS.
ps.gate_polyline_coincidence = "PASS";
```

The plan agrees this is a tautology of CGAL `corefine`, but a
test that *would have caught* a hypothetical CGAL bug or a
property-map miscount is missing entirely. The smoke test
(`test_corefine_smoke.cpp`) does NOT verify, for any pair of
meshes, that the constrained edges of A and B reference the
same 3-D point pairs.

If a future CGAL release changes the constrained-mark
propagation in a way that breaks this invariant, this tool
would silently hand a non-conformal pair to HXT.

**Trigger:** Hypothetical future CGAL bug; or a misuse of
`corefine` (e.g., passing different ECMs by mistake) that
violates the invariant.

**Suggested fix:** Replace the unconditional "PASS" with a real
check that, for every constrained edge in A, an edge with the
same two 3-D endpoints (within `1e-9 m` tolerance — exact
predicates so equality is exact) exists in B's constrained-edge
set:

```diff
@@ tools/corefine_faults/main.cpp lines 339–342
-            // Polyline-edge-coincidence gate is a tautology of corefine
-            // (the constrained edges are by construction the same on
-            // both sides of the cut), so we record PASS.
-            ps.gate_polyline_coincidence = "PASS";
+            // Polyline-edge-coincidence gate: verify, for every
+            // constrained edge in A, that an edge with the same two
+            // 3-D endpoints exists as a constrained edge in B.  This
+            // is documented by CGAL as a tautology of corefine, but a
+            // direct check protects against future CGAL regressions
+            // and against accidental ECM mis-wiring.
+            ps.gate_polyline_coincidence =
+                polyline_edge_coincidence_gate(
+                    meshes[i], ecms[i], meshes[j], ecms[j])
+                ? "PASS" : "FAIL";
+            if (ps.gate_polyline_coincidence != "PASS") {
+                std::cerr << "[corefine_faults] polyline edge-"
+                             "coincidence gate FAILED on pair ("
+                          << ps.short_a << ", " << ps.short_b
+                          << "): a constrained edge in A has no "
+                             "matching constrained edge in B.\n";
+                return 2;
+            }
```

with the helper:
```cpp
bool polyline_edge_coincidence_gate(
    const Mesh& A, const ECMap& ecmA,
    const Mesh& B, const ECMap& ecmB) {
    // Collect canonicalised endpoint pairs (sorted lex) from each side.
    auto collect = [](const Mesh& M, const ECMap& ecm) {
        std::set<std::pair<std::array<double,3>, std::array<double,3>>> s;
        for (auto e : M.edges()) {
            if (!get(ecm, e)) continue;
            auto h = M.halfedge(e);
            const auto& p = M.point(M.source(h));
            const auto& q = M.point(M.target(h));
            std::array<double,3> a{p.x(), p.y(), p.z()};
            std::array<double,3> b{q.x(), q.y(), q.z()};
            if (b < a) std::swap(a, b);
            s.emplace(a, b);
        }
        return s;
    };
    return collect(A, ecmA) == collect(B, ecmB);
}
```

**Test case:**
Add to `test_corefine_smoke.cpp`:
```cpp
void test_polyline_coincidence() {
  Mesh A = make_square(2, 1.0);
  Mesh B = make_square(1, 1.0);
  sc::EdgeConstrainedMap ea, eb;
  install_ecm(A, ea, "e:test");
  install_ecm(B, eb, "e:test");
  sc::corefine_pair(A, B, ea, eb, "A", "B");
  CHECK(polyline_edge_coincidence_gate(A, ea, B, eb),
        "two-square corefine must produce identical 3-D "
        "constrained-edge endpoints on both sides.");
}
```

---

### [R-005] [MODERATE] [tools/corefine_faults/io.cpp:208–235] — `triangle_to_fault.json` writes per-fault `range` from a fresh cumulative cursor, ignoring `per_fault[i].n_triangles` if the in-memory data is inconsistent

**Category:** ASSUMPTION

**Description:**
`write_triangle_to_fault_json` builds the JSON like this:

```cpp
std::int64_t cursor = 0;
for (std::size_t i = 0; i < per_fault.size(); ++i) {
    const auto& f = per_fault[i];
    const std::int64_t lo = cursor;
    const std::int64_t hi = cursor + f.n_triangles;
    cursor = hi;
    os << "    \"" << json_escape(f.short_name) << "\": {"
       << "\"n_triangles\": " << f.n_triangles
       << ", \"range\": [" << lo << ", " << hi << "]}";
    ...
}
os << "  \"n_total_triangles\": " << cursor << '\n';
```

`n_total_triangles` is computed by accumulating `f.n_triangles`
over the loop. The `range` arrays are computed from the same
cumulative cursor. So for any input `per_fault` vector, the
output JSON satisfies `n_total_triangles == sum(faults[*].n_triangles)`
and contiguous ranges by construction.

**This is what makes R-001's "schema validation" tautological**:
the writer cannot produce a JSON that violates the invariants
the validator pretends to check. So even if R-001 is replaced
with a real reread, that check is also tautological unless we
inject a deliberate corruption (e.g., run a separate consistency
check between the in-memory `per_fault` and an INDEPENDENT count
of post-corefine triangles).

**Suggested fix:** Add a *cross-check* between the JSON-loop's
cumulative `cursor` and an INDEPENDENTLY-computed total from
`meshes[*].number_of_faces()`:

```diff
@@ tools/corefine_faults/main.cpp, after writing the JSONs
+    // R-005 cross-check: the total in triangle_to_fault.json must
+    // equal the sum of meshes[i].number_of_faces() — this catches a
+    // mismatch between the writer's loop and the actual mesh state
+    // (e.g., if per_fault was populated from a stale snapshot).
+    {
+        std::int64_t mesh_total = 0;
+        for (const auto& m : meshes)
+            mesh_total += static_cast<std::int64_t>(m.number_of_faces());
+        std::int64_t pf_total = 0;
+        for (const auto& f : per_fault) pf_total += f.n_triangles;
+        if (mesh_total != pf_total) {
+            std::cerr << "[corefine_faults] schema cross-check: "
+                         "in-memory mesh total " << mesh_total
+                      << " != per_fault sum " << pf_total << '\n';
+            return 3;
+        }
+    }
```

This is the actual corruption guard the plan asked for.

**Test case:** Inject inconsistency by mutating one mesh AFTER
populating `per_fault[i].n_triangles`; expect exit 3.

---

### [R-006] [MODERATE] [tools/corefine_faults/corefine.cpp:78–131] — `find_coincident_triangle` returns enumeration counter, not face index; misleading on a mesh with removed faces

**Category:** ASSUMPTION (fragile under mesh mutation)

**Description:**
`find_coincident_triangle` returns `{a_idx, b_idx}` where `a_idx`
/ `b_idx` are counters incremented inside the loop:

```cpp
std::int64_t a_idx = 0;
for (auto f : A.faces()) {
    ...
    table[k].push_back(Entry{a_idx, sorted});
    ++a_idx;
}
```

For a freshly-loaded mesh from `read_ascii_stl` (no removals),
`a_idx` happens to equal `f.idx()` because face indices are
assigned 0, 1, 2, … contiguously. But after any
`remove_face` / `garbage_collect` cycle, this would diverge —
silently. The error message printed in `main.cpp:291` ("fault X
tri 27 …") would point to the wrong triangle.

CGAL's `Surface_mesh::Face_index` exposes `idx()` to get the
underlying integer index, which is stable across removals.

**Suggested fix:**

```diff
@@ tools/corefine_faults/corefine.cpp:88
-    std::int64_t a_idx = 0;
     for (auto f : A.faces()) {
         auto pts = face_points(A, f);
         ...
-        table[k].push_back(Entry{a_idx, sorted});
-        ++a_idx;
+        table[k].push_back(Entry{static_cast<std::int64_t>(f.idx()), sorted});
     }
@@ tools/corefine_faults/corefine.cpp:104
-    std::int64_t b_idx = 0;
     for (auto f : B.faces()) {
         ...
-        if (it != table.end()) {
-            for (const auto& e : it->second) {
-                ...
-                if (same) {
-                    return {e.face_idx, b_idx};
-                }
-            }
-        }
-        ++b_idx;
+        if (it != table.end()) {
+            for (const auto& e : it->second) {
+                ...
+                if (same) {
+                    return {e.face_idx,
+                            static_cast<std::int64_t>(f.idx())};
+                }
+            }
+        }
     }
```

**Test case:**
```cpp
TEST(R006) {
  // Build A with 3 faces, then remove face 1 (collect garbage),
  // then add a NEW face that's coincident with B's face 0.
  // The new face's .idx() will be 3 (not 1), but a counter
  // would report 1 (since the for-loop runs 2 times).
  ...
  auto [ai, bi] = find_coincident_triangle(A, B);
  EXPECT_EQ(ai, 3);   // .idx() of the re-added face
}
```

---

### [R-007] [MODERATE] [tools/corefine_faults/io.cpp:285] — JSON `stl_path` is absolute; not portable across moved output directories

**Category:** QUALITY

**Description:**
The end-to-end run produced:
```json
"stl_path": "/tmp/safs_review_e2e/out/safs_sbmt_millcreek.stl"
```

Paths are stored absolute because `args.out_stl_dir` was
absolute. Moving the output directory breaks consumers (humans
diffing the JSON, comparison harness, archival).

The Python pipeline uses paths relative to the project root
(e.g. `output/two_crossing_2000m/stl_conformal/...`).

**Suggested fix:** Store the path relative to `out_stl_dir` (the
directory the JSON is written to), so moving the directory
preserves correctness:

```diff
@@ tools/corefine_faults/io.cpp:281 inside write_intersection_report_json
-           << ", \"stl_path\": \"" << json_escape(f.stl_path.string()) << "\"}";
+           << ", \"stl_path\": \""
+           << json_escape(
+                  std::filesystem::path(f.stl_path).filename().string())
+           << "\"}";
```

(Just the filename — the JSON lives next to the STL files, so a
filename-only path is unambiguous.)

Update `main.cpp:374` to pass only the filename if the schema
is changed:
```diff
-        per_fault[i].stl_path    = out;
+        per_fault[i].stl_path    = out.filename();
```

---

### [R-008] [LOW] [tools/corefine_faults/corefine.cpp:38] — Missing `<cstring>` include for `std::memcpy`

**Category:** QUALITY

**Description:**
`make_tri_key` calls `std::memcpy` (line 38), but `corefine.cpp`
only includes `<cstdint>`. The build succeeds because
`<unordered_map>` (line 9) transitively pulls in `<cstring>`,
but this is fragile — a future libstdc++ refactor that drops
the transitive include would break the build.

**Suggested fix:**

```diff
@@ tools/corefine_faults/corefine.cpp:6
 #include <array>
 #include <cstdint>
+#include <cstring>
 #include <stdexcept>
```

**Test case:** N/A (compile-only check).

---

### [R-009] [LOW] [tools/corefine_faults/corefine.cpp:24–62] — `make_tri_key` does not normalise +0.0 vs −0.0; potential hash false-negative on coincident triangles

**Category:** EDGE_CASE

**Description:**
`make_tri_key` uses `std::memcpy` to bit-cast doubles to
uint64_t for hashing. IEEE-754 distinguishes `+0.0` from `-0.0`
at the bit level (`0x0000000000000000` vs `0x8000000000000000`),
but `+0.0 == -0.0` when compared as doubles.

If fault A has a vertex at `+0.0` on some axis and fault B has
the same vertex at `-0.0` (e.g., from `ts_to_stl.py`'s clearance
clamp where `(z - clearance) * 0` could yield −0), the hash
mismatch hides the coincident triangle from the preflight.

The probability is low on real CFM data (vertex coordinates are
typically large positive numbers in metres) but non-zero, and
the failure mode is silent: corefine is then called with
coincident triangles, where its behaviour is undefined per CGAL
docs.

**Suggested fix:**

```diff
@@ tools/corefine_faults/corefine.cpp:35
         for (double v : k.sorted_xyz) {
+            if (v == 0.0) v = 0.0;   // canonicalise +0.0 / -0.0
             std::uint64_t bits;
             static_assert(sizeof(bits) == sizeof(v));
             std::memcpy(&bits, &v, sizeof(bits));
             h ^= bits;
             h *= 1099511628211ULL;
         }
```

Apply the same normalisation in the `Point_3` lex sort comparator
in `find_coincident_triangle` (line 95–98 and 109–113).

**Test case:**
```cpp
TEST(R009) {
  Mesh A; A.add_vertex(P3(+0.0, 0.0, 0.0));
          A.add_vertex(P3(1.0, 0.0, 0.0));
          A.add_vertex(P3(0.0, 1.0, 0.0));
          A.add_face(...);
  Mesh B; B.add_vertex(P3(-0.0, 0.0, 0.0));   // bit-different
          B.add_vertex(P3(1.0, 0.0, 0.0));
          B.add_vertex(P3(0.0, 1.0, 0.0));
          B.add_face(...);
  auto [ai, bi] = find_coincident_triangle(A, B);
  CHECK(ai >= 0 && bi >= 0,
        "coincident triangles must match across +0.0/-0.0 jitter");
}
```

---

### [R-010] [LOW] [tools/corefine_faults/io.cpp:68–71] — `repair_polygon_soup(require_same_orientation=false)` could merge front+back faces of a closed surface

**Category:** EDGE_CASE

**Description:**
The implementation passes:
```cpp
PMP::repair_polygon_soup(
    points, polygons,
    CGAL::parameters::erase_all_duplicates(true)
                     .require_same_orientation(false));
```

`require_same_orientation(false)` allows merging two polygons
that share the same vertex set but with opposite winding. For
fault surfaces (open patches with a single defined normal
direction), this is unlikely to fire. But for an STL that
accidentally contains both sides of a closed surface (e.g., a
buggy `ts_to_stl.py` run), CGAL would silently merge them and
produce a non-orientable mesh.

**Suggested fix:** Default to `true` (require same orientation
for duplicate merging); only set to `false` if a future fault
surface format demands it:

```diff
@@ tools/corefine_faults/io.cpp:68–71
     PMP::repair_polygon_soup(
         points, polygons,
         CGAL::parameters::erase_all_duplicates(true)
-                         .require_same_orientation(false));
+                         .require_same_orientation(true));
```

**Test case:** Build a tiny STL with two triangles sharing
3 vertices but with opposite winding; verify with the new flag
that the read fails (or reports the issue) rather than silently
collapsing.

---

### [R-011] [LOW] [tools/corefine_faults/corefine.hpp:59–63] — `corefine_pair` takes ECMs by value (plan said by const reference)

**Category:** DEVIATION

**Description:**
The header signature:
```cpp
CorefineResult corefine_pair(
    Mesh& A, Mesh& B,
    EdgeConstrainedMap ecm_A,                // SHARED across pairs (P-002)
    EdgeConstrainedMap ecm_B,                // SHARED across pairs (P-002)
    std::string_view short_a, std::string_view short_b);
```

The plan (after P-002 was applied) specified `const
EdgeConstrainedMap&`. The implementation takes by value.

This is **functionally fine** — `Property_map` is a CGAL handle
type with value semantics; copying the handle copies a small
struct, not the underlying property storage. But the deviation
makes the "shared across pairs" comment misleading: a reader
might believe the function takes a reference and that mutating
operations propagate. They do, but only because Property_map's
underlying storage is mesh-owned, not because the parameter is a
reference.

**Suggested fix:** Match the plan, eliminating any ambiguity:

```diff
@@ tools/corefine_faults/corefine.hpp:59–63
 CorefineResult corefine_pair(
     Mesh& A, Mesh& B,
-    EdgeConstrainedMap ecm_A,                // SHARED across pairs (P-002)
-    EdgeConstrainedMap ecm_B,                // SHARED across pairs (P-002)
+    const EdgeConstrainedMap& ecm_A,         // SHARED across pairs (P-002)
+    const EdgeConstrainedMap& ecm_B,         // SHARED across pairs (P-002)
     std::string_view short_a, std::string_view short_b);
```

And matching update in `corefine.cpp:142–146`. CGAL's
`PMP::corefine` accepts the parameter as a property-map alias,
so passing `ecm_A` (a `const&`) works fine.

---

## Summary

- Critical issues: **2** (R-001 tautological schema check, R-002
  plan acceptance not met)
- Moderate issues: **5** (R-003 silent gate failure, R-004
  unverified polyline-coincidence, R-005 missing cross-check,
  R-006 enumeration vs face_index, R-007 absolute STL path)
- Low issues: **4** (R-008 missing include, R-009 +0/-0 hash,
  R-010 repair flag, R-011 by-value ECM)
- Plan compliance: **PARTIAL**. Phase 0 acceptance is met
  end-to-end (build + hello + smoke test pass). Phase 1
  acceptance is *partially* met: the binary runs, all four
  smoke tests pass, but (a) R-002 fails the stated 500-floor,
  (b) R-001's schema check doesn't actually validate, (c) R-003
  silently swallows manifold-gate failures.
- Verdict: **PASS WITH FIXES**. The build infrastructure is
  sound and the recipe is correct. R-001 + R-003 are
  silent-failure paths that must be closed before Phase 2 lands
  on top of them; R-002 needs a plan-vs-data reconciliation
  even if the answer is "the plan was over-counting".

## Suggested fix order

1. **R-001 + R-005** (real schema check + cross-check) — both
   change the same area; one fix.
2. **R-003** (manifold-gate exit propagation) — three-line edit;
   immediately closes a silent failure.
3. **R-002** (plan-vs-actual reconciliation) — investigate
   first: measure the polyline length from the CGAL output and
   compare to the Python output. If they agree, update the plan
   to 143 ± 10. If they disagree, file a follow-up to debug
   CGAL's intersection emit.
4. **R-004** (real polyline-coincidence check) — adds a real
   gate where the current code asserts a tautology.
5. **R-006 + R-007** (face_index correctness + relative STL
   paths) — independent quality fixes.
6. **R-008, R-009, R-010, R-011** — cleanups; can land
   together.

## Unreviewed Areas

- **`PMP::corefine` `throw_on_self_intersection(true)` behaviour
  on actually-self-intersecting input.** The smoke test does not
  exercise this path; manual verification was deferred. Add a
  smoke test that constructs a self-intersecting fault and
  asserts the runtime_error fires (caught by main.cpp:329 with
  exit code 2).
- **End-to-end pass through `generate_safs_mesh.py` with HXT
  on the CGAL output.** Phase 1 only asserts the conformal STLs
  are produced; Phase 3 will verify HXT accepts them. If HXT
  rejects, R-002's "143 vs 500" question becomes moot because
  the topology is wrong regardless of count.
- **Behaviour on the disjoint-faults case** (no triangle pairs
  intersect): not exercised by the smoke test. Plan §Phase 1
  edge cases requires a no-op with empty `pairs[…]` entry.
- **CGAL header shim correctness** (`tools/cgal_patch/`). I
  read the README and CMakeLists wiring; the shim's `this->g
  != nullptr` substitution is documented as semantically
  equivalent. I did not diff the shim file against the
  upstream conda-forge header to confirm it's a 3-line patch
  vs a wholesale rewrite. The provenance section of the README
  asserts the diff is exactly three sed substitutions; trust
  but verify before the next CGAL bump.
