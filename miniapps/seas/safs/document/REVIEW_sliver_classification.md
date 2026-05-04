# Code Review: Sliver classification (2026-04-30)

## Review Scope
- Plan: none formal — review prompted by user observation that the
  earlier classification ("144/167 single-fault curvature wedges,
  0 cross-fault wedges") contradicted the visible mesh: slivers are
  visible at every fault-fault intersection in the .vtu, regardless
  of dihedral angle.
- Files reviewed: `/tmp/classify_unfixable.py` (the diagnostic script
  that produced the original analysis), and the chain of conclusions
  it drove ("Option 2a rejected", "γ_min floor is structural",
  "single-fault curvature is the root cause").
- Domain context: `mesh/output/clean6_2000m_cgal/output/fault_provenance.json`
  schema (key `triangle_indices_in_msh` is indexed against the FULL
  .msh element list, not the tag-100 subset); user's screenshot of
  visible slivers along fault-fault intersection lines.

## Findings

### [R-001] CRITICAL [classify_unfixable.py:main] — Index-space mismatch silently zeroes the per-fault component map

**Category:** BUG

**Description:**
The script builds `fault_tri_list` as the 0-based subset of tag-100
triangles (range 0..5397). It then reads `triangle_indices_in_msh`
from `fault_provenance.json` and stores them into `tri_to_fault`
indexed against the **same** 0..5397 space, gated by

```python
for idx in idxs:
    if 0 <= idx < len(fault_tri_list):
        tri_to_fault[idx] = k
```

But `triangle_indices_in_msh` is indexed against the **full .msh
element list** — verified ranges per fault are
`[9724, 10148]`, `[10149, 10655]`, ..., `[13264, 15121]`. Every
provenance index ≥ 9724, every `len(fault_tri_list) == 5398`, so
the gate `idx < 5398` is **false for every provenance index**. None
are written; `tri_to_fault` stays the default `-1` everywhere.

**Trigger:**
Running `classify_unfixable.py` against any mesh whose provenance
JSON uses full-msh indexing (which is the actual schema).

**Actual behavior:**
- `tri_to_fault[i] == -1` for every fault triangle.
- For every sliver, `comps = {-1}` regardless of which faults its
  vertices touch.
- `len(comps) <= 1` is **always** true → every multi-component sliver
  is mis-classified as single-component.
- Output: `single-comp = 144, multi-comp = 0`. **All 141
  cross-fault wedge slivers are hidden.**

**Expected behavior:**
After fix: `single-comp = 3, multi-comp = 141`. The dominant
sliver mode is cross-fault wedges at five specific fault-pair
intersections (sbmt_saf × sbmt_missioncreek = 41, sbmt_saf ×
sbmt_millcreek = 39, mjvs_saf × sbmt_saf = 27, coav_missioncreek ×
mult_ssaf_banning = 15, sbmt_millcreek × sbmt_missioncreek = 13).

**Suggested fix:**
Convert provenance indices into the 0-based tag-100 subset by
subtracting the offset of the smallest provenance index (which is
where the tag-100 block starts in the .msh element ordering).

```diff
 if prov_path.exists():
     prov = json.loads(prov_path.read_text())
     faults_dict = prov.get('faults', {})
     if faults_dict:
+        all_idxs = []
+        for entry in faults_dict.values():
+            all_idxs.extend(entry.get('triangle_indices_in_msh', []))
+        offset = min(all_idxs) if all_idxs else 0
         tri_to_fault = np.full(len(fault_tri_list), -1, dtype=int)
         for k, (name, entry) in enumerate(faults_dict.items()):
             names.append(name)
             idxs = entry.get('triangle_indices_in_msh', [])
             for idx in idxs:
-                if 0 <= idx < len(fault_tri_list):
-                    tri_to_fault[idx] = k
+                j = idx - offset
+                if 0 <= j < len(fault_tri_list):
+                    tri_to_fault[j] = k
```

**Test case:**
```python
def test_R001_provenance_index_offset_resolves_multi_component():
    # Run classify_unfixable.py against the clean6_2000m_cgal mesh.
    # Before fix: multi-comp == 0.  After fix: multi-comp >= 100.
    out = run_classify("mesh/output/clean6_2000m_cgal/output/safs_clean6_2000m.msh")
    assert out["single_comp"] == 3
    assert out["multi_comp"] == 141
```

---

### [R-002] CRITICAL [classify_unfixable.py:main] — No validation that the provenance map actually populated anything

**Category:** ASSUMPTION

**Description:**
After parsing provenance, the script never checks how many entries
of `tri_to_fault` got assigned a non-`-1` value. With R-001, all
5398 entries stay `-1` and the script proceeds silently to use
`tri_to_fault[tri_idx]` as a fault id. A single sanity assertion
would have caught the bug at first run.

**Trigger:**
Any time the provenance JSON schema differs from what the script
assumes (schema drift, version skew, off-by-one).

**Actual behavior:**
Silently produces classification dominated by single-component
slivers (because every `comps` set evaluates to `{-1}`).

**Expected behavior:**
Hard-fail (or warn loudly) if `tri_to_fault` has zero non-`-1`
entries — a classification with 100% unknown-fault triangles is
useless and points at a schema mismatch.

**Suggested fix:**
Add a self-check immediately after building `tri_to_fault`:

```diff
             for k, (name, entry) in enumerate(faults_dict.items()):
                 ...
                 for idx in idxs:
                     j = idx - offset
                     if 0 <= j < len(fault_tri_list):
                         tri_to_fault[j] = k
+            n_assigned = int(np.sum(tri_to_fault >= 0))
+            if n_assigned == 0:
+                raise RuntimeError(
+                    f"provenance index→fault map is empty: "
+                    f"{len(all_idxs)} provenance indices, "
+                    f"{len(fault_tri_list)} tag-100 triangles, "
+                    f"none in range. Schema mismatch?")
+            if n_assigned < len(fault_tri_list):
+                print(f"WARN: {len(fault_tri_list) - n_assigned} "
+                      f"fault triangles have no provenance entry "
+                      f"(of {len(fault_tri_list)})", file=sys.stderr)
```

**Test case:**
```python
def test_R002_empty_provenance_map_raises():
    # Synthetic provenance with all-out-of-range indices.
    fake_prov = {"faults": {"f0": {"triangle_indices_in_msh":
                                   [99999, 100000]}}}
    # Expect RuntimeError, not silent zero-fill.
    with pytest.raises(RuntimeError, match="provenance index"):
        run_classify_with_prov(fake_prov)
```

---

### [R-003] CRITICAL [conversation analysis] — Root-cause conclusion was based on R-001's bad output

**Category:** DEVIATION (consequence of R-001)

**Description:**
The earlier debug session's conclusion — "144/167 unfixable
slivers are single-fault curvature wedges; the floor is structural"
— was driven entirely by R-001's mis-classification. The
corrected data shows **141/144 four-fault slivers are
cross-fault polyline wedges**, concentrated at five specific
fault-pair intersections (the SAF strands meeting Mill Creek /
Mission Creek branches). Cross-fault wedges are **not** an
"inherent property of conformal Delaunay against a 2-manifold";
they are a property of the cross-fault polyline DIHEDRAL ANGLES
emerging from autorefine.

This invalidates two downstream conclusions:
1. The recommendation "accept γ_min = 1.36×10⁻⁴ as a structural floor"
   was wrong — it's not structural.
2. The negative-result interpretation of Option 2a (`Field[2].SizeMin
   = res_f * 0.3`) — that "the wedge tet's 4 vertices are all on
   the same fault, so bulk refinement cannot reach it" — was
   geometrically wrong. The wedges span multiple faults; bulk
   refinement IN the dihedral cavity between them MIGHT in principle
   help, though the empirical bit-identity of γ_min between baseline
   and Option 2a still indicates the specific Field[2].SizeMin knob
   was not the right one. Different geometric mechanisms (2b, 3, 4)
   need re-evaluation against the corrected diagnosis.

**Trigger:**
Trusting the buggy classification output and propagating it through
multiple Option-evaluation rounds without re-checking.

**Actual behavior:**
Five rounds of analysis built on a false premise; user finally
caught it from visual inspection of the .vtu.

**Expected behavior:**
The classification should have been validated against an independent
signal (e.g., visual confirmation that intersection-line slivers
match the multi-comp count) before being used to reject mechanisms.

**Suggested fix (process, not code):**

1. Re-run `classify_unfixable.py` after fixing R-001 (already done
   in this session — output above).
2. Update the `safs.geo` comment block (lines 124–134) that documents
   the rejected Option 2a — its rationale ("the wedge tet's 4
   vertices are all on the same fault") is now known to be wrong.
   The empirical observation (γ_min unchanged at 1.358e-04 with
   `SizeMin = res_f * 0.3`) is still valid, but the GEOMETRIC
   EXPLANATION for why must be replaced — likely "the dihedral
   cavity at fault-fault intersections is too narrow for bulk
   vertices spaced at 600 m".
3. Re-evaluate Options 2b, 3, 4 (in REVIEW the user's earlier table)
   against the corrected diagnosis — cross-fault dihedral wedges
   are amenable to local Steiner-point insertion (CGAL `exude_mesh_3`
   / TetGen `-rq`) along the polyline, which is a much more
   surgical operation than what we discussed under the "single-fault
   curvature" framing.
4. The `code-debug` log written above the safs.geo comment
   ("Sliver-fix attempt rejected (2026-04-30...)") should not be
   left as a permanent record without correcting the geometric
   reason it gives.

**Test case:**
```python
def test_R003_corrected_classification_is_persisted():
    # The downstream artifacts (safs.geo comment block, future plan
    # docs) must reflect the corrected counts.  Failing this test
    # means the bad analysis is still load-bearing somewhere.
    geo = open("mesh/safs.geo").read()
    # The old (wrong) rationale must be replaced or annotated.
    assert "4 vertices on the same fault" not in geo, (
        "safs.geo still cites the R-001-driven incorrect rationale")
```

---

### [R-004] MODERATE [classify_unfixable.py:main] — `prov_path` derivation is brittle

**Category:** EDGE_CASE

**Description:**
`prov_path = msh_path.parent / "fault_provenance.json"` assumes
the provenance JSON sits next to the mesh — true for the
`output/clean6_*_cgal/output/` layout but not for ad-hoc runs.
When missing, the script silently falls back to "no per-fault
classification" without any complaint, so a user who points the
tool at a mesh in a different directory gets the same broken
single-comp result as R-001 — but for a different reason and
without any signal.

**Trigger:**
Running the script with a mesh whose provenance JSON is in a
non-default location, or against a mesh produced by a pipeline
without provenance.

**Actual behavior:**
`tri_to_fault` is `None` → control flow skips the multi-comp
classification entirely → single-comp/multi-comp counters report
0/0 for the 4-fault group.

**Expected behavior:**
Either accept a `--provenance-json` CLI argument, or warn loudly
when the default path doesn't exist.

**Suggested fix:**
```diff
 def main():
     msh_path = Path(sys.argv[1])
-    prov_path = msh_path.parent / "fault_provenance.json"
+    prov_path = msh_path.parent / "fault_provenance.json"
+    if not prov_path.exists():
+        print(f"WARN: no provenance JSON at {prov_path}; "
+              f"per-fault classification will be skipped",
+              file=sys.stderr)
```

**Test case:**
```python
def test_R004_missing_provenance_warns():
    # Move provenance aside; expect a stderr WARN, not silent zero.
    pass
```

---

### [R-005] LOW [classify_unfixable.py:main] — Single-component count includes one-vertex-on-fault cases

**Category:** QUALITY

**Description:**
The classifier records `comps.add(int(tri_to_fault[tri_idx]))`
for each fault triangle a vertex belongs to. A vertex on a CROSS-
FAULT POLYLINE belongs to triangles from BOTH faults, so its
`comps` set contains TWO ids. This is the path we rely on for
multi-comp detection — and it works after R-001 is fixed. But the
counting also means that a sliver where vertex v is the ONLY
multi-fault vertex (the other 3 vertices belong to one fault each)
is classified as multi-comp via that single vertex. That is the
correct geometric meaning, but the report's "cross-fault wedge"
language could mislead readers into expecting all 4 vertices to
straddle two faults.

**Suggested fix:**
Clarify the report wording — distinguish "tet has at least one
polyline vertex" from "tet straddles two faults". Cosmetic only.

**Test case:** none required (LOW).

---

## Summary
- Critical issues: 3 (R-001, R-002, R-003)
- Moderate issues: 1 (R-004)
- Low issues: 1 (R-005)
- Plan compliance: N/A (no formal plan)
- Verdict: **FAIL** — R-001 is a correctness bug that invalidated
  every downstream conclusion in the prior debug session. Must fix
  R-001 + R-002 + R-003 before any further sliver analysis.

## Unreviewed Areas
- The `safs.geo` comment block at lines 124–134 that documents the
  rejected Option 2a uses the now-incorrect "4 vertices on the same
  fault" rationale. R-003 calls this out for fix; the comment needs
  to be edited but only after the geometric reason is properly
  re-derived.
- The "Option 2 / 2a / 2b / 3 / 4 / 5" decision matrix in the prior
  conversation table needs to be re-scored against the corrected
  diagnosis (141 cross-fault wedges, not 144 single-fault curvatures).
  Multi-comp wedges may be tractable by local Steiner insertion in
  the dihedral cavity, which was previously deemed unviable under
  the wrong diagnosis. Out of scope for THIS review (this review is
  about the classifier bug, not the next sliver-fix plan).
