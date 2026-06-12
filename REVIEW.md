# Code Review: 2026-06-12 — PLAN_mesh_quality_safv4_remesh_2026-06-12.md (plan-document review)

> Previous REVIEW.md content (free-surface slice, Round 4, 2026-06-02) is superseded;
> git history retains it (commit f959238).

## Review Scope
- Plan under review: `miniapps/seas/safs/project_7.0_preferred/document/PLAN_mesh_quality_safv4_remesh_2026-06-12.md`
  (the plan itself is the artifact; findings are spec defects that /code-fix should apply AS EDITS TO THE PLAN
  before /code-implement runs).
- Evidence consulted: `meshing/code/run_vtu_quality_model_SAFv4_2km_topo.log`, `meshing/code/run_quality.log`,
  `document/STATUS_tetgen_pipeline.md`, `document/PLAN_revert_to_gmsh.md`, archive commit `577646e`
  (`git diff-tree`), archived tool CLIs (`git show 577646e:...`), and NEW measurements run during this review
  on the per-surface VTUs in `meshing/vtu/` (fault-boundary contact, ribbon planarity, shell rim coincidence).
- Domain context: repo CLAUDE.md, `miniapps/seas/CLAUDE.md` (Gmsh v2.2 limitation), session memory
  (`safv4-gocad-mesh-diagnosis-remesh-plan`).

## Findings

### [R-001] [CRITICAL] [PLAN Phase 2/3:fault border handling] — Faults touch the BOTTOM surface and EACH OTHER; plan only handles the DEM contact

**Category:** ASSUMPTION / DEVIATION (from physical geometry)

**Description:**
The plan's extend-and-clip treatment exists ONLY for the fault top border against the DEM
(`extend_fault_above_dem.py`, Phase 4 clip). Phase 3 explicitly models the other borders as free
("Fault tip inside the domain ... allowed; the tip border is a free border"). Measured this review
(coordinate coincidence at 1 mm on `meshing/vtu/*.vtu`):

```
SAF (MJVS)  minus: 53 nodes coincident with `bottom` (fault z_min = -19,329 = bottom plane)
Garnet      minus:  5 nodes coincident with `bottom`
MJVS x MULT  (SAF x Banning)      :  9 coincident nodes
MJVS x SBMT  (SAF x Garnet Hill)  : 21 coincident nodes
MULT x SBMT  (Banning x Garnet)   : 68 coincident nodes
(no fault touches any ribbon)
```

SAF and Garnet Hill are trimmed AT the domain bottom, and all three faults are trimmed against each
other at junction lines (sealed GOCAD model). After Phase 3 remeshing with free borders, those
borders drift off the contacted surface by up to the remesh displacement: some vertices end up
epsilon-OUTSIDE the bottom plane or epsilon-through the host fault (poke), others epsilon-short
(gap). Poke-through past the bottom plane = fault facet outside the closed PLC shell -> tetgen PLC
error (the plan's own highest-listed risk). Junction gap/poke = fragmented corefine polylines +
splinter triangles violating the 100 m floor and the max_diff_m = 0 conformality gate.

**Trigger:**
Any Phase 3 run on the real SAFv4 geometry (the contacts exist today; see measurements).

**Actual behavior (per current plan):**
Bottom-trimmed and junction-trimmed borders are remeshed free; nothing re-establishes contact;
Phase 4's clip never removes below-bottom or through-host-fault overhangs.

**Expected behavior:**
Every fault border edge is classified by contact target and handled so the final fault patch
terminates EXACTLY on the contacted entity (or is a genuine interior tip).

**Suggested fix (plan edits):**
1. Phase 2, Detailed Requirements item 6 — extend the manifest:
```diff
- and for each fault the **trace polyline**: ordered border-edge chains whose both endpoints are
- (post-weld) shared with the DEM node set.
+ and for each fault a **border classification**: every border edge assigned to exactly one of
+ {dem (trace), bottom, ribbon_<name>, fault_<basename> (junction), interior_tip}, by post-weld
+ shared-node test against each candidate surface (both endpoints shared -> contact). Emit ordered
+ chains per class. Acceptance must reproduce the measured contacts: SAF-bottom 53 nodes,
+ Garnet-bottom 5, MJVS-MULT 9, MJVS-SBMT 21, MULT-SBMT 68 (1 mm tol).
```
2. Phase 3 — generalize the extension (rename `extend_fault_above_dem.py` ->
   `extend_fault_borders.py`):
```diff
- For each trace chain: extrude each trace vertex by +z h_ext ...
+ For each border chain, extrude by class:
+   dem (trace):    +z by h_ext (default 1500 m), verified strictly above the DEM (unchanged).
+   bottom:         -z by h_ext_bot (default 500 m) below z_bottom = -19329.4 m (bottom is exactly
+                   planar, max plane-fit deviation 0.000 m, measured).
+   fault junction: along the local fault-plane outward direction (average of border-edge in-plane
+                   normals) by h_ext_junc (default 750 m), so the extended patch crosses the host
+                   fault transversally; the host fault is NOT extended at the shared junction.
+                   Host/guest assignment: the fault whose border lies ON the other's interior is
+                   the guest (from the Phase 2 classification).
+   interior_tip:   no extension.
+   ribbon:         none observed in SAFv4 (assert; abort if Phase 2 classification finds one).
```
3. Phase 3 corefine pairs note: junction pairs (MJVS-MULT, MJVS-SBMT, MULT-SBMT) are now transversal
   crossings; the intersection polyline replaces the old trimmed junction (same pattern as the trace).
4. Phase 4 clip: see R-002 (the clip predicate must remove ALL overhangs, not only above-DEM).

**Test case:**
```python
def test_R001_border_classification_reproduces_measured_contacts():
    # 1mm quantized key sets per surface from meshing/vtu/*.vtu (script run 2026-06-12, output above):
    # assert len(saf_keys & bottom_keys) == 53
    # assert len(garnet_keys & bottom_keys) == 5
    # assert len(mjvs_keys & mult_keys) == 9
    # assert len(mjvs_keys & sbmt_keys) == 21
    # assert len(mult_keys & sbmt_keys) == 68
    # Phase 2 manifest border classes must be non-empty for exactly these pairs + dem.
```

---

### [R-002] [CRITICAL] [PLAN Phase 4:clip_fault_at_dem.py] — Clip predicate is DEM-heightfield-only; bottom/junction overhangs survive into the tetgen input

**Category:** BUG (spec-level)

**Description:**
Phase 4 item 3 drops fault triangles with `centroid z > z_DEM(x,y) + 1e-6` — only the above-DEM
overhang. With R-001's required extensions (and even without them, with remesh border drift), fault
triangles can lie BELOW the bottom plane or PAST a host fault's surface at a junction. A
below-bottom triangle has centroid z < z_DEM trivially, so the DEM test KEEPS it; it lies outside
the closed boundary shell, and tetgen rejects the PLC (or meshes a phantom region).

**Trigger:**
A fault patch extended/drifted past z = -19,329.4 m (SAF, Garnet — both bottom-trimmed today), or
past a host fault at a junction (MULT-SBMT has 68 contact nodes today).

**Actual behavior:** out-of-domain / through-host-fault triangles reach `tetgen_mesh.py`.

**Expected behavior:** the tetgen input contains fault triangles only strictly inside the closed
boundary shell, terminating on the corefined intersection polylines.

**Suggested fix (plan edit, Phase 4 item 3, rename script `clip_fault_overhangs.py`):**
```diff
- **Clip** (`clip_fault_at_dem.py`): drop every fault-marked triangle whose centroid is above
- the DEM: above-test = centroid z > z_DEM(x, y) + 1e-6 ...
+ **Clip** (`clip_fault_overhangs.py`), applied per fault patch AFTER autorefine:
+   (a) shell containment: drop every fault-marked triangle whose centroid is OUTSIDE the closed
+       boundary shell. Test = ray parity (+z ray from centroid, count crossings against
+       boundary-marked triangles; outside if even). The DEM-heightfield shortcut is valid for the
+       top only and misses below-bottom overhangs — do not use it as the sole predicate.
+   (b) junction overhang: for each junction pair from the Phase 2 classification, flood-fill the
+       extended (guest) fault's triangles from its extension strip across shared edges, stopping
+       at the corefined intersection polyline; drop the filled (overhang) component.
+   (c) keep the 1e-6 margin + any-vertex-on-surface tie-break rule for triangles incident to the
+       intersection polylines (unchanged).
```
Acceptance additions: `0 fault triangles outside the shell (re-run test (a) as a gate)`;
`every junction border edge lies on its corefined junction polyline`.

**Test case:**
```python
def test_R002_clip_drops_below_bottom_triangle():
    # Synthetic: closed unit shell (flat top z=0 as 'dem', bottom z=-1, 4 sides);
    # one vertical fault quad spanning z in [-1.2, +0.2] (extended both ways).
    # Old predicate (z > z_dem) drops only the z>0 part -> triangle at z=-1.1 KEPT  -> FAIL.
    # New predicate (a) drops both overhangs -> remaining fault strictly inside [-1, 0] -> PASS.
```

---

### [R-003] [MODERATE] [PLAN Numerical constraints / Phase 5 gates] — Hard gate G2a (0 tets with eta <= 0.05) contradicts the only empirical baseline of the chosen pipeline

**Category:** ASSUMPTION

**Description:**
The plan's own evidence base (`STATUS_tetgen_pipeline.md`) reports the proven pipeline left
~70 tets with q < 0.05 (0.21%) on the 2000 m flat-box fixture, and documents that every attempted
mitigation was marginal (tetgen optim: 70 -> 60) with MMG3D banned. The plan asserts the root cause
(non-uniform fault density 100 m - 2.5 km) disappears with uniform 500 m faults — plausible but
unproven. As written, G2a is a HARD Phase 5 gate, so if the assertion is even slightly wrong, every
production run mechanically fails into Phase 6 with no calibrated expectation, on a 5-20 M-tet run
that costs hours per iteration.

**Trigger:** first 3-fault coarse rung (c) producing any eta <= 0.05 tet.

**Actual behavior (per plan):** G2a binary-fails the ladder; Phase 6 entered with no agreed budget.

**Expected behavior:** the gate is calibrated where iteration is cheap (rung (c)) and the
hard/soft split is an explicit user decision BEFORE the production run.

**Suggested fix (plan edit, Phase 5 run ladder item 4c + gates):**
```diff
   c. **3-fault coarse**: full SAFv4 at 2000 m fault target.
+     Rung (c) is also the G2 calibration point: record count(eta<=0.05) and count(eta<=0.1) with
+     the Phase-1 localization report. If count(eta<=0.05) > 0, STOP and present the report with a
+     recommendation (retune Phase-3 sizing vs accept-localized) BEFORE rung (d); G2a applies to
+     rung (d) as agreed at that decision point. G2a remains the default target.
```

**Test case:**
```python
def test_R003_rung_c_emits_calibration_record():
    # After rung (c): gate JSON must contain keys
    # {"n_eta_le_0p05", "n_eta_le_0p1", "sliver_z_hist", "sliver_boundary_attribution"}
    # and the ladder must refuse to start rung (d) when n_eta_le_0p05 > 0
    # without an explicit --accept-g2a-waiver flag.
```

---

### [R-004] [MODERATE] [PLAN Phase 1:fault_embedding] — Specified algorithm (global face->tet incidence map) is infeasible at 31 M tets in Python

**Category:** BUG (spec-level / feasibility)

**Description:**
`fault_embedding` is specified as "builds a face->tet incidence map (sorted-triple key over all
4*M tet faces, chunked)". For M = 31.3e6 that is ~125e6 dict entries; a Python dict of tuple keys
costs ~100+ B/entry -> >12 GB, exceeding the measured ~8 GB envelope of the existing tools.
Chunking the tet loop does not bound the MAP size — the map is global. The tool must run on the
known-bad 31 M-tet `.inp` (Phase 1 acceptance) and on the final production mesh.

**Trigger:** Phase 1 acceptance run on `model_SAFv4_mesh_2km_topo.inp`.

**Actual behavior:** memory blow-up / swap death on the baseline mesh.

**Expected behavior:** memory bounded by the fault size (~1e4-1e5 triangles), not the tet count.

**Suggested fix (plan edit, Phase 1 item 1):**
```diff
- def fault_embedding(model: SurfaceModel, fault_name: str) -> dict
-     # builds a face->tet incidence map (sorted-triple key over all 4*M tet faces, chunked);
+ def fault_embedding(model: SurfaceModel, fault_name: str) -> dict
+     # builds the key set of FAULT triangles only (sorted node triples, ~1e4-1e5 entries),
+     # then streams tets in 2e6 chunks: form each chunk's 4 face triples (vectorized np.sort on
+     # the (4*chunk, 3) index array), pack each triple into one int64/typed key, and count
+     # matches against the fault key set (np.isin); accumulate per-fault-triangle incidence.
+     # Memory O(n_fault_tri + chunk), never O(n_tet).
```

**Test case:**
```python
def test_R004_embedding_memory_bound():
    # synthetic 2e6-tet mesh with a 1e3-triangle fault: peak traced allocation of
    # fault_embedding < 1.5x the mesh-array footprint (tracemalloc), and incidence == 2
    # for all interior fault triangles.
```

---

### [R-005] [MODERATE] [PLAN Phase 0 + Phase 4 --box regression] — Regression fixture inputs were EXCLUDED from the archive; Phase 0 restore pathspec cannot support the mandated regression

**Category:** DEVIATION / ASSUMPTION

**Description:**
Phase 4 acceptance requires re-running "the archived 2000 m flat-box fixture" and reproducing
"the same triangle count 8,202". But commit `577646e`'s message explicitly EXCLUDES the corefined
STLs ("data_corefined/*.stl ... regenerable from raw_data/*.ts"), and Phase 0 restores ONLY
`code_preprocess/` — not `raw_data/` (which IS in the commit). After Phase 0 as written: no `.ts`
inputs, no corefined STLs, fixture unregenerable. Additionally, regenerating via
`corefine_faults.py --res 2000` under a possibly-drifted CGAL 6.1.1 env is not guaranteed
bit-identical, so "same triangle count 8,202" is an unverifiable equality.

**Trigger:** Phase 4 acceptance run (and Phase 5 rung 4a).

**Actual behavior:** regression step fails for missing inputs; even when regenerated, exact-count
equality may fail spuriously.

**Expected behavior:** restore everything the regression needs; compare structurally.

**Suggested fix (plan edits):**
```diff
 Phase 0:
-   `git checkout 577646e -- miniapps/seas/safs/project_7.0_preferred/code_preprocess`
+   `git checkout 577646e -- miniapps/seas/safs/project_7.0_preferred/code_preprocess \
+        miniapps/seas/safs/project_7.0_preferred/raw_data \
+        miniapps/seas/safs/project_7.0_preferred/data_corefined`   # manifests; STLs regenerated
+   Add step: regenerate the 2000 m fixture (`python code_preprocess/corefine_faults.py --res 2000`)
+   and record its manifest as the regression reference; diff against the archived May manifest and
+   report any drift (CGAL env changes show up here, before any new code).
 Phase 4 acceptance:
- [ ] `--box` regression: ... reproduces the STATUS_tetgen_pipeline.md soup (same triangle count
-     8,202 +- the 8 documented dedup losses downstream; markers identical).
+ [ ] `--box` regression (structural): merged soup self-intersection-free; 6 fault markers present;
+     per-marker triangle counts within 2% of the regenerated fixture's corefine output; downstream
+     tetgen rung (Phase 5 4a): embedding 100%, coverage 100/100, n_tets within 10% of 33,650.
```

**Test case:**
```bash
# after Phase 0 (current plan): demonstrates the gap
ls miniapps/seas/safs/project_7.0_preferred/raw_data/*.ts            # FAILS (not restored)
git diff-tree --name-only -r 577646e | grep -c 'raw_data'            # > 0 (exists in the commit)
```

---

### [R-006] [MODERATE] [PLAN Phase 2 item 6 / acceptance] — Weld-count acceptance number 2,939 is scoped to FAULT nodes only; the global weld will exceed it

**Category:** BUG (spec-level)

**Description:**
The baseline "2,939 coincident groups" was measured over the 6,490 fault-surface nodes ONLY.
Phase 2's weld is global (all 5.16 M nodes), and this review additionally measured
boundary-boundary coincidences (NE30-ribbon x DEM: 142 nodes; ribbons x bottom: 24-69 each;
ribbon-ribbon corner chains: 3-4) — so the global group count is strictly larger than 2,939. The
acceptance "n_groups_merged = expected 2,939-vs-baseline check" fails on CORRECT output (or worse,
gets "fixed" by restricting the weld to fault nodes, breaking shell watertightness).

**Trigger:** first Phase 2 run.

**Actual behavior:** spurious acceptance failure on correct output.

**Expected behavior:** compare like with like.

**Suggested fix (plan edit, Phase 2 item 6):**
```diff
- Also record the weld stats (n_groups_merged = expected 2,939-vs-baseline check)
+ Also record the weld stats, split by class: n_groups_total (informational; expected > 2,939
+ because boundary rims are also duplicated/coincident), and n_groups_fault_subset (groups
+ containing >= 1 fault-surface node) which MUST equal the baseline 2,939 at tol = 1e-3 m.
```

**Test case:**
```python
def test_R006_weld_stats_scoping():
    # run extract_safv4_surfaces.py on the frozen .inp;
    # assert manifest["weld"]["n_groups_fault_subset"] == 2939
    # assert manifest["weld"]["n_groups_total"] > 2939   # boundary rims exist (measured)
```

---

### [R-007] [MODERATE] [PLAN Phase 5 item 2 (.msh tag map)] — Garbled, self-contradicting block inside the normative tag-map spec

**Category:** QUALITY (actively misleading in a normative section)

**Description:**
The `.msh` output spec contains an abandoned half-written map followed by an inline editorial aside:

```
 rock volume                  = 1
 fault_SAF (MJVS)             = 101   (alphabetical basename order: Banning, Garnet, SAF
 fault_Banning                = 102    -> Banning=101? NO — fix the map explicitly:)
```

before the actual normative map. The parenthetical alphabetizes by SHORT name (Banning, Garnet,
SAF) while the Constraints section mandates alphabetical FULL-basename order — which the normative
map does satisfy (SAFS-SAFZ-MJVS < SAFS-SAFZ-MULT < SAFS-SAFZ-SBMT). A literal implement agent can
read the artifact as licence to assign Banning=101. One spec, stated once.

**Trigger:** /code-implement reading Phase 5 item 2.

**Actual behavior:** two conflicting maps + an editing artifact in a normative block.

**Expected behavior:** single normative map with the ordering rule stated.

**Suggested fix (plan edit):**
```diff
-     rock volume                  = 1
-     fault_SAF (MJVS)             = 101   (alphabetical basename order: Banning, Garnet, SAF
-     fault_Banning                = 102    -> Banning=101? NO — fix the map explicitly:)
-     ```
-     Normative default tag map (alphabetical, stable):
+     Normative default tag map (faults ordered alphabetically by FULL CFM basename:
+     "SAFS-SAFZ-MJVS-..." < "SAFS-SAFZ-MULT-..." < "SAFS-SAFZ-SBMT-..."):
```
(keep the normative block that follows unchanged).

**Test case:**
```bash
grep -c 'NO — fix the map explicitly' \
  miniapps/seas/safs/project_7.0_preferred/document/PLAN_mesh_quality_safv4_remesh_2026-06-12.md
# must be 0 after fix
```

---

### [R-008] [MODERATE] [PLAN Phase 3/Phase 4 boundary] — Intersection-resolution responsibility is ambiguous; Phase 4 text invites deferring fault x DEM to autorefine, which has NO min-edge cleanup

**Category:** ASSUMPTION / DEVIATION risk

**Description:**
Phase 3 owns corefine + polyline cleanup (the only machinery that achieved >= 96 m edges and
max_diff_m = 0). Phase 4 then lists "extension strip x DEM" as an example of "residual
intersections" for autorefine to resolve. The strip x DEM intersection IS the trace — if it reaches
Phase 4 unresolved, autorefine splits it but nothing collapses the resulting sub-100 m splinter
triangles BELOW the trace; the clip removes only the above-DEM side, so G1 fails late, at the most
expensive stage, with no attribution.

**Trigger:** implementer routes fault x DEM through autorefine (the Phase 4 example invites it),
or a Phase 3 cleanup gap goes undetected.

**Actual behavior:** sub-100 m splinters can survive into the tetgen input.

**Expected behavior:** Phase 3 MUST fully resolve and clean all fault x fault and fault x DEM
(incl. strips) intersections; Phase 4 autorefine is a safety net that must find nothing new
involving a fault marker.

**Suggested fix (plan edit, Phase 4 item 2):**
```diff
- `PMP::autorefine_triangle_soup` resolves residual intersections (e.g. extension strip x DEM,
-  strip x ribbon) with marker propagation ...
+ `PMP::autorefine_triangle_soup` is a SAFETY NET: all fault x fault and fault x DEM (incl.
+  extension strips) intersections are already corefined + cleaned in Phase 3. Expected residual
+  work here: strip x ribbon clamps and exact-duplicate removal only. Gate: autorefine's
+  resolved-intersection count involving any FAULT marker must be 0; if > 0, STOP and report which
+  Phase 3 pair leaked (a Phase 3 bug — not something to absorb here).
```

**Test case:**
```python
def test_R008_autorefine_fault_intersections_zero():
    # run Phase 4 on Phase 3 output; parse autorefine_merged --verbose log;
    # assert n_resolved_intersections_with_fault_marker == 0
    # negative control: feed an UNcorefined fault+dem pair -> count > 0 and exit nonzero.
```

---

### [R-009] [MODERATE] [PLAN Phase 3:Files to Modify] — The graded DEM remesh has no implementation home; "modify corefine_set.cpp only if needed" understates a new C++ deliverable

**Category:** DEVIATION (missing file-level plan item)

**Description:**
The archived isotropic remeshing is uniform (one `--mesh-edge-size` per run);
`distance_sizing_field.h` targets CGAL Mesh_3 (`mesh_volume.cpp`), not PMP::isotropic_remeshing.
The DEM step needs: (a) a sizing field h(d) over distance-to-trace driving
`PMP::isotropic_remeshing` (sizing-field overload, CGAL >= 5.6, present in 6.1), (b) rim border
edges constrained with `protect_constraints = true` AND rim vertices excluded from relocation
(`vertex_is_constrained_map`), or the watertight rim weld breaks silently. Neither exists in any
restored file; an implement agent reaches Phase 3 and discovers a missing tool mid-phase.

**Trigger:** Phase 3 DEM step.

**Actual behavior:** no tool to run; improvisation risk inside the most delicate phase.

**Expected behavior:** explicit deliverable with interface.

**Suggested fix (plan edit, Phase 3 Files to Create):**
```diff
 ### Files to Create
 - `meshing/code/extend_fault_above_dem.py`
+ - `code_preprocess/corefine_cgal/remesh_graded.cpp` — graded isotropic remeshing of one surface:
+   `remesh_graded IN.stl OUT.stl --trace-polylines manifest.json --h-near 500 --h-far 2500
+    --d-near 1000 --d-far 9000 [--protect-border] [--verbose]`
+   Implements h(d) = h_near + (h_far - h_near) * clamp((d - d_near)/(d_far - d_near), 0, 1) as a
+   model of PMPSizingField, with distance d to the trace polyline vertices via a CGAL KD/AABB
+   query; border edges constrained (edge_is_constrained_map, protect_constraints = true) and
+   border vertices fixed (vertex_is_constrained_map) — the rim vertex set must be byte-identical
+   pre/post. Add the target to corefine_cgal/CMakeLists.txt.
```

**Test case:**
```python
def test_R009_rim_immutability_and_grading():
    # remesh a synthetic 20km x 20km bumpy grid with a straight 'trace' down the middle:
    # (1) border vertex coordinate set identical pre/post (sorted byte compare);
    # (2) median edge in [400, 700] m for faces with centroid within 1 km of the trace;
    # (3) median edge > 1500 m for centroids beyond 9 km.
```

---

### [R-010] [LOW] [PLAN Phase 1 acceptance] — Circular dependency: Phase 1 acceptance references a Phase 5 artifact

**Category:** QUALITY

**Description:** "check_mesh_quality.py on a `.msh` produces the same bulk metrics as on the
equivalent `.vtu` ... (e.g., a Phase 5 smoke output)" — but Phase 1 must be accepted before
Phase 5 exists (Phase 1 is "Required by: Phases 2, 3, 4, 5").

**Suggested fix:**
```diff
- (cross-check on any small fixture, e.g. a Phase 5 smoke output)
+ (cross-check on an existing v2.2 fixture already in the repo, e.g. a tpv205/tpv102 `.msh`
+  under miniapps/seas/, converted to `.vtu` via meshio)
```

---

### [R-011] [LOW] [PLAN Phase 5 run ladder] — Coarse rungs reuse production-hardcoded sizing constants

**Category:** EDGE_CASE

**Description:** DEM grading (500/2500 m, 1/9 km window), `--polyline-spacing 250`, and tetgen
`--lc-near 500` are stated only for the 500 m production target. Rungs (b)/(c) at 2000 m need all
of these scaled, or the coarse rungs do not exercise the code paths they are meant to de-risk.

**Suggested fix:**
```diff
+ All resolution-dependent knobs are functions of the rung target R: fault mesh_edge_size = R,
+ polyline_spacing = R/2, DEM h_near = R, h_far = 5R, tetgen lc-near = R, lc-far = 10R,
+ G1f bands scaled by R/500. Record resolved values in each rung's gate JSON.
```

---

### [R-012] [LOW] [PLAN gate G3c] — No operational definition of "fault-trace node" on the final `.msh`

**Category:** EDGE_CASE

**Description:** G3c says "every fault-trace node IS a DEM-surface node (same node ID)", but the
final `.msh` carries no trace marker; the checker needs a rule for finding trace nodes.

**Suggested fix:**
```diff
+ Operational G3c on .msh: compute fault border edges (exactly 1 incident fault-tagged triangle,
+ per fault tag). Classify each: shared with a top(201)-tagged triangle edge -> trace; with
+ bottom(202)/sides(203) -> termination; shared with another fault tag -> junction; else tip.
+ G3c passes iff every border edge is classified (no orphans) and trace-edge endpoints appear in
+ top-surface triangles with the SAME node IDs. Report per-class counts per fault.
```

---

## Summary
- Critical issues: 2 (R-001, R-002 — both stem from the unverified "faults only touch the DEM"
  assumption, falsified by measurement this review: SAF/Garnet are bottom-trimmed and all three
  faults touch each other at junction lines)
- Moderate issues: 7 (R-003 .. R-009)
- Low issues: 3 (R-010 .. R-012)
- Plan compliance: n/a (the plan is the artifact under review; against the user's three goals it
  is directionally sound and evidence-based, but Phases 2-4 are mis-specified for the real fault
  border topology)
- Verdict: **FAIL — must fix before proceeding to /code-implement.** R-001/R-002 change the
  Phase 2 manifest schema, the Phase 3 extension tool, and the Phase 4 clip; implementing the plan
  as written would hit a tetgen PLC failure (or silent non-conformal junctions) at the most
  expensive rung of the ladder.

## Unreviewed Areas
- The archived C++ sources (`corefine_set.cpp`, `autorefine_merged.cpp`, `polyline_cleanup.h`, ...)
  were not line-audited; the review trusts the STATUS doc's empirical record. Phase 0 builds them
  before any new code depends on them.
- CGAL 6.1 API availability for the PMP sizing-field overload and constrained-vertex maps was
  asserted from documentation knowledge, not compiled — R-009's test compiles it at Phase 0/3.
- meshio `gmsh22` physical-tag fidelity — covered only by the plan's existing MFEM read smoke gate.
- z-extents of the fault-fault junction chains (whether junctions reach the DEM) were not mapped;
  the R-001 Phase 2 classification produces this as output.
