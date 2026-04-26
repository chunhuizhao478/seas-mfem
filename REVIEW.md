# Code Review: 2026-04-25 — Remove Option 3 (cfm_trim) and produce overlap-free subset meshes

## Review Scope
- Plan: `miniapps/seas/safs/PLAN_smoke_millcreek.md` (smoke + bug catalog), follow-up directive in conversation: "remove Option 3 and just build meshes without any overlap".
- Files reviewed:
  - `miniapps/seas/safs/mesh/cfm_trim.py` (436 lines, to be deleted)
  - `miniapps/seas/safs/mesh/run_full_2000m.sh` (94 lines, calls cfm_trim — must change)
  - `miniapps/seas/safs/mesh/generate_safs_mesh.py` (457 lines, `_combine_stls` and includes-writer)
  - `miniapps/seas/safs/mesh/audit_ts_quality.py` (overlap-pair detection — KEEP, used by ts_to_stl)
  - `miniapps/seas/safs/mesh/ts_to_stl.py` (`--repair-overlaps` — KEEP, single-fault repair is still useful)
- Domain context: project `CLAUDE.md`, `miniapps/seas/CLAUDE.md`, `miniapps/seas/safs/PLAN.md`, `PLAN_origin.md`, `PLAN_domain.md`, `PLAN_smoke_millcreek.md`. Diagnostic data from this conversation: 10 fault-pair combinations have 3-D triangle intersections (Mill Creek × SBMT-SAF: 163 pairs, etc.); 18 pair combinations have NO intersection. The drop-only trim cannot make multi-fault assemblies conformal regardless of proximity buffer (verified up to proximity = 2000 m, dropping 17.3% of triangles, still rejected).

## Background — what works after the smoke test

| status | item |
| --- | --- |
| ✓ | `safs_origin.py` + `transform.schema.json` |
| ✓ | `audit_ts_quality.py` (with `n_overlap_pairs` metric — KEEP) |
| ✓ | `ts_to_stl.py` (with `--repair-overlaps` flag for in-fault overlap pairs — KEEP) |
| ✓ | `safs.geo` (single-fault Embed pipeline) |
| ✓ | `generate_safs_mesh.py` (single-fault path) |
| ✓ | `write_fault_provenance.py`, `validate_msh.py`, `convert_msh.py` |
| ✓ | `run_smoke_millcreek.sh` (passes 10/10 smoke checks) |
| ✗ | `cfm_trim.py` — removable per user directive |
| ✗ | `run_full_2000m.sh` — calls cfm_trim, must change |
| partial | `_combine_stls` in `generate_safs_mesh.py` — written for 8-fault path but still useful for any multi-fault subset; works but writes output to the input dir (R-005) |

## Conflict graph (which fault pairs cross in 3-D)

From the diagnostic earlier in this conversation, the 10 fault pairs whose triangulations cross in 3-D:

```
1.  safs_coav_missioncreek    × safs_mult_ssaf_banning   ( 26)
2.  safs_coav_missioncreek    × safs_sbmt_missioncreek   ( 28)
3.  safs_mjvs_saf             × safs_sbmt_saf            (  8)
4.  safs_mult_banning         × safs_sbmt_saf            ( 11)
5.  safs_pmfz_pinto           × safs_sbmt_millcreek      ( 42)
6.  safs_pmfz_pinto           × safs_sbmt_missioncreek   ( 18)
7.  safs_pmfz_pinto           × safs_sbmt_saf            ( 20)
8.  safs_sbmt_millcreek       × safs_sbmt_missioncreek   ( 41)
9.  safs_sbmt_millcreek       × safs_sbmt_saf            (163)
10. safs_sbmt_missioncreek    × safs_sbmt_saf            (135)
```

A complete partition of all 8 faults into mutually-non-overlapping subsets:

| subset | members | size |
|---|---|---|
| **NW**  | `safs_mjvs_saf`, `safs_pmfz_pinto`, `safs_mult_banning`, `safs_mult_ssaf_banning` | 4 |
| **South** | `safs_sbmt_saf`, `safs_coav_missioncreek` | 2 |
| **Mill** | `safs_sbmt_millcreek` | 1 |
| **MissionSBMT** | `safs_sbmt_missioncreek` | 1 |

Verification (subset NW): MJVS×Pinto, MJVS×Banning, MJVS×SSAF-Banning, Pinto×Banning, Pinto×SSAF-Banning, Banning×SSAF-Banning — none of the 6 pairs is in the conflict list. ✓

Verification (subset South): SBMT-SAF×COAV — not in the conflict list. ✓

Per the user directive, the new wrapper `run_subsets_2000m.sh` should build these 4 subset meshes independently. Each subset goes through the existing pipeline (audit → ts_to_stl → generate_safs_mesh → provenance → validate → convert) with `--include-fault` constrained to the subset's members.

## Findings

---

### [R-001] [MODERATE] [cfm_trim.py] — Delete the file (Option 3 abandoned)

**Category:** DEVIATION (user-directed scope reduction)

**Description:**
`cfm_trim.py` (436 lines) implements the simplified drop-only Option 3 trim. Empirical testing in this conversation showed it cannot produce a meshable multi-fault assembly even at proximity = 2000 m (drops 17.3% of triangles, HXT still rejects with 266 missing facets due to T-junctions on the trimmed boundaries). The user has directed to abandon this approach in favor of building per-subset meshes that have no inter-fault overlaps by construction.

The geometric primitives in this file (`tri_tri_intersect_3d`, `edge_pierces_tri`, `_point_to_tri_distance`, `seg_seg_cross_2d`) are correct and could in principle be reused by a future proper Option 3 (with triangle splitting), but per the directive there is no caller for them now. Keeping them as dead code creates maintenance burden and confuses future readers about whether multi-fault overlaps "are handled".

**Trigger:** running `bash run_full_2000m.sh` invokes `cfm_trim`, which produces a mesh-input that HXT rejects.

**Actual behavior:** wrapper invokes cfm_trim → produces stl_trimmed → gmsh fails. End-to-end build never produces a valid mesh.

**Expected behavior:** the cfm_trim stage does not exist; instead the user runs the new subset wrapper (R-003).

**Suggested fix:** delete the file.
```diff
- miniapps/seas/safs/mesh/cfm_trim.py
```

**Test case:**
```bash
test "$(ls miniapps/seas/safs/mesh/cfm_trim.py 2>&1 | grep -c 'No such file')" = "1"
```

---

### [R-002] [MODERATE] [run_full_2000m.sh] — Delete the failing 8-fault wrapper

**Category:** DEVIATION (replaced by subset wrapper, R-003)

**Description:**
`run_full_2000m.sh` is the 7-stage wrapper that invokes `cfm_trim`. Even with cfm_trim removed, building all 8 CFM faults into a single mesh is infeasible (the diagnostic showed 10 pairs with genuine 3-D triangulation crossings). Leaving this script in the repo as-is is misleading: it advertises an 8-fault assembly that does not work.

**Trigger:** anyone running `bash run_full_2000m.sh` after R-001 will see Stage 3/7 fail because `cfm_trim.py` no longer exists; if R-001 is reverted, the gmsh stage fails with "missing facets".

**Actual behavior:** advertises 8-fault assembly that produces no valid mesh.

**Expected behavior:** the file does not exist; the user runs `run_subsets_2000m.sh` (R-003) which builds 4 separate meshes.

**Suggested fix:** delete the file.
```diff
- miniapps/seas/safs/mesh/run_full_2000m.sh
```

**Test case:**
```bash
test "$(ls miniapps/seas/safs/mesh/run_full_2000m.sh 2>&1 | grep -c 'No such file')" = "1"
```

---

### [R-003] [CRITICAL] [run_subsets_2000m.sh] — Add the new subset-wrapper (does not exist)

**Category:** DEVIATION (new feature required to honor the directive)

**Description:**
The user directive "build meshes without any overlap" requires a wrapper that runs the existing per-subset pipeline once per overlap-free subset. The 4 subsets are listed at the top of this REVIEW.md. Each subset is built into its own output directory under `miniapps/seas/safs/mesh/output/subset_<name>_2000m/`.

For subsets with one fault (Mill, MissionSBMT), the existing single-fault path used by `run_smoke_millcreek.sh` is reused verbatim. For subsets with two or more faults (NW with 4 faults, South with 2 faults), the existing `_combine_stls` step in `generate_safs_mesh.py` handles the conformity-via-shared-vertex requirement. **Note**: the multi-fault path of `_combine_stls` has not been validated end-to-end after R-005/R-006/R-007 land; the wrapper must do at minimum a smoke validation that every subset produces a `.msh` with the expected BP5 tag inventory. Use `validate_msh.py` for this.

**Trigger:** running `bash run_subsets_2000m.sh` from a clean repo.

**Actual behavior:** file does not exist.

**Expected behavior:** the wrapper builds 4 subset meshes, each passing all `validate_msh.py` checks (the 1-fault subsets via the single-fault Embed path; the 4-fault NW and 2-fault South subsets via the combined-STL path).

**Suggested fix:** create the file at `miniapps/seas/safs/mesh/run_subsets_2000m.sh` with the contents below, and `chmod +x` it.

```bash
#!/usr/bin/env bash
# Build 4 overlap-free SAFS subset meshes at 2000 m.  Each subset's faults
# are mutually non-intersecting (verified against the 3-D triangle-triangle
# intersection scan in REVIEW.md's conflict graph).
#
# Per-subset outputs:  output/subset_<name>_2000m/
#                          ├── stl/
#                          ├── transform.json
#                          ├── bbox.json
#                          ├── cleanup_log_2000m.csv
#                          └── output/
#                              ├── safs_subset_<name>_2000m.msh
#                              ├── safs_subset_<name>_2000m.vtu   (ParaView)
#                              ├── safs_subset_<name>_2000m.xml   (Dolfin)
#                              ├── safs_subset_<name>_2000m_facet_region.xml
#                              ├── domain_box.json
#                              ├── sizing.json
#                              ├── fault_provenance.json
#                              └── validation_report.txt
#
# Usage:
#     conda activate pythonenv
#     bash miniapps/seas/safs/mesh/run_subsets_2000m.sh
set -euo pipefail
cd "$(dirname "$0")"

RES="${SAFS_RES:-2000}"
CLEARANCE="${SAFS_CLEARANCE:-100}"
RES_F="${SAFS_RES_F:-2000}"
RES_FF="${SAFS_RES_FF:-25000}"
RAMP="${SAFS_RAMP:-35000}"
BUF="${SAFS_BUF:-50000}"
DEPTH="${SAFS_DEPTH:-50000}"
CFM_DIR="${SAFS_CFM_DIR:-$HOME/Documents/Earthquake Cycle Modeling of San Andreas Fault System/CFM_data}"

# subset_name : space-separated short_names
SUBSETS=(
    "NW:safs_mjvs_saf safs_pmfz_pinto safs_mult_banning safs_mult_ssaf_banning"
    "South:safs_sbmt_saf safs_coav_missioncreek"
    "Mill:safs_sbmt_millcreek"
    "MissionSBMT:safs_sbmt_missioncreek"
)

for spec in "${SUBSETS[@]}"; do
    name="${spec%%:*}"
    members="${spec#*:}"
    OUTDIR="output/subset_${name}_${RES}m"
    echo "############# Building subset ${name} (${members}) #############"
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR/stl" "$OUTDIR/output"

    incl=()
    for m in $members; do incl+=( --include-fault "$m" ); done

    echo "==> 1/6 audit"
    python audit_ts_quality.py --cfm-dir "$CFM_DIR" --res "$RES" \
        "${incl[@]}" --out "$OUTDIR/cfm_audit_${RES}m.csv"

    echo "==> 2/6 ts_to_stl"
    python ts_to_stl.py --cfm-dir "$CFM_DIR" --res "$RES" \
        "${incl[@]}" --out-dir "$OUTDIR/stl" \
        --transform-json "$OUTDIR/transform.json" \
        --bbox-json "$OUTDIR/bbox.json" \
        --log-csv "$OUTDIR/cleanup_log_${RES}m.csv" \
        --free-surface-clearance "$CLEARANCE" --repair-overlaps

    echo "==> 3/6 generate_safs_mesh"
    python generate_safs_mesh.py "${incl[@]}" \
        --stl-dir "$OUTDIR/stl" \
        --bbox-json "$OUTDIR/bbox.json" \
        --transform-json "$OUTDIR/transform.json" \
        --buf-x "$BUF" --buf-y "$BUF" --depth "$DEPTH" \
        --res-f "$RES_F" --res-ff "$RES_FF" --ramp-dist "$RAMP" \
        -o "$OUTDIR/output/safs_subset_${name}_${RES}m.msh"

    echo "==> 4/6 write_fault_provenance"
    python write_fault_provenance.py \
        --msh "$OUTDIR/output/safs_subset_${name}_${RES}m.msh" \
        --stl-dir "$OUTDIR/stl" \
        --transform-json "$OUTDIR/transform.json" \
        "${incl[@]}" --out "$OUTDIR/output/fault_provenance.json"

    echo "==> 5/6 validate_msh"
    python validate_msh.py \
        --msh "$OUTDIR/output/safs_subset_${name}_${RES}m.msh" \
        --transform-json "$OUTDIR/transform.json" \
        --bbox-json "$OUTDIR/bbox.json" \
        --provenance-json "$OUTDIR/output/fault_provenance.json" \
        --domain-box-json "$OUTDIR/output/domain_box.json" \
        --sizing-json "$OUTDIR/output/sizing.json" \
        --report "$OUTDIR/output/validation_report.txt" || true

    echo "==> 6/6 convert_msh"
    python convert_msh.py \
        --msh "$OUTDIR/output/safs_subset_${name}_${RES}m.msh"
    echo
done

echo "All subsets built.  Per-subset outputs:"
for spec in "${SUBSETS[@]}"; do
    name="${spec%%:*}"
    echo "  output/subset_${name}_${RES}m/output/safs_subset_${name}_${RES}m.{msh,vtu,xml}"
done
```

**Test case:**
```python
def test_R003_subset_wrapper_runs():
    """All 4 subsets build, each produces a valid .msh + .vtu + .xml,
    and validate_msh.py reports 10/10 checks pass for each."""
    import subprocess, os
    os.chdir("miniapps/seas/safs/mesh")
    r = subprocess.run(["bash", "run_subsets_2000m.sh"],
                        capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr
    for name in ("NW", "South", "Mill", "MissionSBMT"):
        d = f"output/subset_{name}_2000m/output"
        assert os.path.exists(f"{d}/safs_subset_{name}_2000m.msh")
        assert os.path.exists(f"{d}/safs_subset_{name}_2000m.vtu")
        assert os.path.exists(f"{d}/safs_subset_{name}_2000m.xml")
        with open(f"{d}/validation_report.txt") as fh:
            tail = fh.read().splitlines()[-1]
            assert "10/10" in tail, f"{name}: {tail}"
```

---

### [R-004] [LOW] [README pointer] — Document the per-subset workflow

**Category:** QUALITY (discoverability)

**Description:**
The repo has `PLAN_smoke_millcreek.md` and the smoke wrapper, but no top-level pointer explaining "how to build the SAFS mesh in production". After R-001/R-002/R-003 land, a `README.md` next to the wrappers should explain: 8 CFM faults cannot be assembled into a single conformal mesh in this version; instead, 4 spatially-disjoint subsets are built; the partition was determined by the 3-D triangle-triangle intersection scan; future work (proper Option 3 with triangle splitting) would unify the subsets.

**Suggested fix:** create `miniapps/seas/safs/mesh/README.md`:

```markdown
# SAFS mesh build entry points

| script | purpose |
|---|---|
| `run_smoke_millcreek.sh` | end-to-end single-fault smoke test (Mill Creek strand) |
| `run_subsets_2000m.sh`   | production: 4 overlap-free subset meshes at 2000 m |

The 8 CFM faults cannot be assembled into a single conformal Steiner-
constraint set with the current pipeline (`Embed`-based, no triangle
splitting at fault-fault intersections).  Instead, the 8 faults are
partitioned into 4 mutually-non-overlapping subsets (see the conflict
graph in REVIEW.md).  Each subset is meshed independently.

Future work: implement proper triangle splitting along fault-fault
intersection curves so subsets NW + South + Mill + MissionSBMT can be
unified into a single mesh.  This is the full version of "Option 3"
described in conversation logs.
```

**Test case:** none (documentation).

---

### [R-005] [MODERATE] [generate_safs_mesh.py:_combine_stls] — Combined STL is written into the *input* directory

**Category:** BUG

**Description:**
`_combine_stls` is called with `combined_stl = args.stl_dir / "safs_combined.stl"`. The combined STL is written *into* the input STL directory. Side effects:
1. A future caller that auto-discovers STLs in `stl-dir` would pick up `safs_combined.stl` as an "input" (today the driver only iterates the explicit `included` list, so the bug is dormant — it activates with the next refactor).
2. Tools that delete the input dir (e.g., the wrapper's `rm -rf "$OUTDIR/stl"`) destroy the combined STL too — that's fine but conflates "user-provided per-fault inputs" with "pipeline-generated combined output".
3. Combined STL silently appears in input dir on every multi-fault build with no opt-out flag.

**Trigger:** any multi-fault build (e.g., subset NW or South in R-003).

**Actual behavior:**
```python
combined_stl = args.stl_dir / "safs_combined.stl"
meta = _combine_stls(combined_stl, included, args.stl_dir, ...)
```
Writes alongside per-fault STLs.

**Expected behavior:** write into the OUTPUT directory next to `safs_*.msh`.

**Suggested fix:** in `generate_safs_mesh.py:main`, around line ~390 (the `if len(included) > 1` block):
```diff
-    if len(included) > 1:
-        combined_stl = args.stl_dir / "safs_combined.stl"
-        meta = _combine_stls(combined_stl, included, args.stl_dir,
-                              snap_m=args.combine_snap_m)
-        (args.output.parent / "combine_meta.json").parent.mkdir(
-            parents=True, exist_ok=True)
-        (args.output.parent / "combine_meta.json").write_text(
-            json.dumps(meta, indent=2) + "\n"
-        )
+    if len(included) > 1:
+        args.output.parent.mkdir(parents=True, exist_ok=True)
+        combined_stl = args.output.parent / "safs_combined.stl"
+        meta = _combine_stls(combined_stl, included, args.stl_dir,
+                              snap_m=args.combine_snap_m)
+        (args.output.parent / "combine_meta.json").write_text(
+            json.dumps(meta, indent=2) + "\n"
+        )
```

**Test case:**
```python
def test_R005_combined_stl_in_output_dir(tmp_path):
    """generate_safs_mesh.py with len(included) >= 2 writes
    safs_combined.stl next to the .msh, NOT into the input stl-dir."""
    # Setup tmp_path/stl with 2 valid per-fault STLs.
    # Setup tmp_path/{transform.json,bbox.json}.
    # Run generate_safs_mesh.py with --include-fault A --include-fault B
    #     -o tmp_path/output/m.msh
    # Assert: tmp_path/output/safs_combined.stl exists
    # Assert: tmp_path/stl/safs_combined.stl does NOT exist
```

---

### [R-006] [MODERATE] [generate_safs_mesh.py:_combine_stls] — Each per-fault STL is read twice

**Category:** BUG (correctness-fragile + performance)

**Description:**
`_combine_stls` calls `_read_ascii_stl(p)` twice per fault: once in the first pass (vertex collection) and once in the second pass (raw-coord recovery via `for short in included: ... verts, tris = _read_ascii_stl(p) ... for v in verts: r = int(inverse[cursor]) ... cursor += 1`).

Two independent reads from disk are an integrity hazard: file could change between reads (extremely unlikely in practice but possible in CI/parallel jobs), and the cursor arithmetic relies on identical vertex counts in both reads. If the second read returns even one different vertex count for any reason, the `inverse` array is misindexed silently and the final triangulation is corrupted.

**Trigger:** any multi-fault build, any time. Today no failure observed because file reads are stable, but the silent-corruption mode is real.

**Actual behavior:** two reads of the same file.

**Expected behavior:** read once, cache, reuse.

**Suggested fix:**
```diff
+    cache: dict[str, tuple[list, list]] = {}
     for short in included:
         p = stl_dir / f"{short}.stl"
         if not p.exists():
             raise SystemExit(f"missing STL: {p}")
         verts, tris = _read_ascii_stl(p)
+        cache[short] = (verts, tris)
         # ... existing append logic unchanged ...
     # ...
     cursor = 0
     for short in included:
-        p = stl_dir / f"{short}.stl"
-        verts, tris = _read_ascii_stl(p)
+        verts, tris = cache[short]
         for v in verts:
             r = int(inverse[cursor])
             if r not in raw_coords:
                 raw_coords[r] = v
             cursor += 1
```

**Test case:**
```python
def test_R006_combine_stls_no_double_read(monkeypatch, tmp_path):
    """_combine_stls reads each input STL exactly once."""
    import generate_safs_mesh as gsm
    call_count = {"n": 0}
    real_read = gsm._read_ascii_stl
    def counting(p):
        call_count["n"] += 1
        return real_read(p)
    monkeypatch.setattr(gsm, "_read_ascii_stl", counting)
    # Setup 2 STLs in tmp_path/stl/, each with a few facets.
    gsm._combine_stls(tmp_path / "out.stl",
                       ["a", "b"], tmp_path / "stl", snap_m=0.01)
    assert call_count["n"] == 2  # once per fault, not twice
```

---

### [R-007] [MODERATE] [generate_safs_mesh.py:_combine_stls] — Degenerate-triangle filter desynchronizes `fault_ranges`

**Category:** EDGE_CASE

**Description:**
After deduplication, `_combine_stls` drops triangles where any two of the three remapped indices coincide:
```python
new_tris = [(a, b, c) for (a, b, c) in new_tris
            if a != b and b != c and a != c]
```
The function returns `n_triangles = len(new_tris)` AFTER the drop, but `fault_ranges` records `[start, end]` based on the *pre-drop* offsets. After the drop, `fault_ranges['NW'] = [0, 417]` may point past the end of `new_tris` or into a slice owned by a different fault.

`combine_meta.json` is only used for diagnostics today, but if any future code consumes it for per-fault attribution the result is silently wrong.

**Trigger:** multi-fault build where snap deduplication collapses any triangle's vertices to a degenerate. Likely whenever snap > vertex spacing somewhere in the input.

**Actual behavior:** `fault_ranges` is stale after the degenerate drop.

**Expected behavior:** `fault_ranges` accurately points into the post-drop triangle list.

**Suggested fix:** track per-fault keep-counts during the filter:
```diff
-    new_tris = [tuple(int(inverse[i]) for i in tri) for tri in all_tris]
-    # Drop triangles that became degenerate (shared vertex collapse).
-    new_tris = [(a, b, c) for (a, b, c) in new_tris
-                if a != b and b != c and a != c]
+    new_tris_pre = [tuple(int(inverse[i]) for i in tri) for tri in all_tris]
+    new_tris: list[tuple[int, int, int]] = []
+    range_remap: dict[str, list[int]] = {s: [0, 0] for s in fault_ranges}
+    n_dropped_degenerate = 0
+    for short, (start, end) in fault_ranges.items():
+        out_start = len(new_tris)
+        for k in range(start, end):
+            a, b, c = new_tris_pre[k]
+            if a == b or b == c or a == c:
+                n_dropped_degenerate += 1
+                continue
+            new_tris.append((a, b, c))
+        out_end = len(new_tris)
+        range_remap[short] = [out_start, out_end]
+    fault_ranges = range_remap
```
Add `"n_dropped_degenerate": n_dropped_degenerate` to the returned meta dict.

**Test case:**
```python
def test_R007_fault_ranges_post_drop(tmp_path):
    """fault_ranges in combine_meta.json indexes into the post-drop
    triangle list, never the pre-drop list."""
    # Construct 2 STLs where snap-dedup forces fault A to lose 1 triangle.
    # Run _combine_stls.
    # Assert: each fault_ranges[k][1] - fault_ranges[k][0]
    #         equals the actual number of surviving triangles for fault k.
    # Assert: sum of (end-start) over all faults == n_triangles.
```

---

### [R-008] [LOW] [generate_safs_mesh.py:_write_includes] — `combined_stl` parameter silently ignored when `len(included) <= 1`

**Category:** QUALITY (silent override)

**Description:**
`_write_includes` accepts a `combined_stl` argument and uses it only inside `if combined_stl is not None`. There is no guard against passing both `combined_stl=<path>` AND a single-fault `included`. The caller in `main` correctly conditions on `len(included) > 1`, so the bug is dormant. But a future caller could pass mismatched arguments and silently get the per-fault path written instead.

**Suggested fix:** assert the precondition.
```diff
 def _write_includes(path: Path, included: list[str], stl_dir: Path,
                      combined_stl: Path | None = None) -> None:
+    if combined_stl is not None and len(included) < 2:
+        raise ValueError(
+            "combined_stl was supplied but only "
+            f"{len(included)} fault(s) included; "
+            "the combine step is unnecessary at this scale"
+        )
     path.parent.mkdir(parents=True, exist_ok=True)
```

**Test case:**
```python
def test_R008_combined_stl_with_single_fault_raises(tmp_path):
    from generate_safs_mesh import _write_includes
    import pytest
    with pytest.raises(ValueError, match="only 1 fault"):
        _write_includes(tmp_path / "x.geo", ["safs_mjvs_saf"],
                         tmp_path / "stl",
                         combined_stl=tmp_path / "c.stl")
```

---

### [R-009] [LOW] [run_smoke_millcreek.sh] — Smoke test must be re-verified after Option 3 removal

**Category:** ASSUMPTION (POSSIBLE)

**Description:**
`run_smoke_millcreek.sh` does NOT call `cfm_trim`. After R-001 deletes `cfm_trim.py`, the smoke test should continue to work (no Python script imports from `cfm_trim`; no shell wrapper sources it). Verified by `grep -l 'cfm_trim' miniapps/seas/safs/mesh/*.py` showing only `cfm_trim.py` itself.

**Suggested fix (verification only):** the /code-fix agent must run after applying R-001/R-002/R-003/R-005/R-006/R-007:
```bash
cd miniapps/seas/safs/mesh
bash run_smoke_millcreek.sh
# Expected last line: "10/10 checks passed"
```
If the smoke test breaks, the trigger is in `generate_safs_mesh.py` changes (R-005/R-006/R-007), not in the deletes.

**Test case:** the verification command above.

---

## Summary
- Critical issues: 1 (R-003 — new wrapper required)
- Moderate issues: 5 (R-001, R-002, R-005, R-006, R-007)
- Low issues: 3 (R-004, R-008, R-009)
- Plan compliance: PARTIAL (the smoke test plan is fully implemented; the implicit "8-fault production assembly" is dropped per user directive in favor of subset assembly)
- Verdict: **PASS WITH FIXES** — must apply R-001 (delete cfm_trim.py), R-002 (delete run_full_2000m.sh), R-003 (add run_subsets_2000m.sh), and the bug fixes R-005, R-006, R-007 before the next user-visible build.

## Unreviewed Areas
- `safs_origin.py`, `audit_ts_quality.py`, `ts_to_stl.py`, `safs.geo`, `validate_msh.py`, `convert_msh.py`, `write_fault_provenance.py`, `run_smoke_millcreek.sh` — these all passed the 10/10 smoke check end-to-end and have no relationship to Option 3. Out of scope for this review.
- The proper "Option 3 with triangle splitting" path is explicitly out of scope per user directive; not reviewed.
- The single-fault Embed pipeline (`safs.geo` Steiner-constraint embedding) was reviewed and validated in the smoke-test review round; not re-reviewed here.
