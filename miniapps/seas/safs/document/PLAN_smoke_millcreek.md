# Implementation Plan: Single-Fault Smoke Test — Mill Creek Strand

## Overview
End-to-end smoke test of the SAFS mesh pipeline (audit → STL → Gmsh → tagged `.msh` →
validate) using exactly one CFM fault — the **Mill Creek fault strand**
(`SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4`, strike 286°, dip 83°, 444 vrtx / 794
trgl at the 2 km resolution). Goal: produce a runnable mesh that BP5 can load, and
surface every bug the multi-fault pipeline is likely to hit *before* we wire all eight
CFMs together. No code outside `miniapps/seas/safs/mesh/` is touched.

This plan is testing-only. It does **not** replace `PLAN.md`, `PLAN_origin.md`, or
`PLAN_domain.md` — it exercises them.

## Mill Creek strand: properties relevant to this test

```
File:         SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4_2000m.ts
N vertices:   444
N triangles:  794
Avg strike:   286°    (NW–SE oriented)
Avg dip:       83°    (near-vertical)
Surface area: 1089.64 km²
UTM bbox:     X ∈ [469 788, 528 501] m
              Y ∈ [3 767 227, 3 784 534] m
              Z ∈ [-18 014,    2 593] m
Local bbox:   X ∈ [-30 212,  +28 501] m   (after origin -500 000 E)
              Y ∈ [+2 227,  +19 534] m   (after origin -3 765 000 N)
              Z ∈ [-18 014, 0] m         (after z>0 → 0 clamp)
z>0 vrtxs:    45  (10.1% — must be flattened in cleanup)
```

After Phase 2 cleanup the surface should still be 2-manifold (this fault's CFM file is
known-clean per the master plan's audit step).

## Constraints

- **Reuse, do not refactor.** The smoke test invokes the same scripts the multi-fault
  pipeline will use. If a script needs a flag to operate on a single fault, add the
  flag — do not fork.
- **Test artefacts go under `mesh/output/smoke_millcreek/`.** Do not pollute the
  default output directory the multi-fault pipeline will write to.
- **Default origin from `PLAN_origin.md`.** UTM (500 000, 3 765 000, 0) — do not pass
  overrides. The smoke test exists *because* this is the production-path origin.
- **One BP5-tag scheme.** All Mill Creek triangles end up under Physical Surface 100.
  No per-fault sub-tags. This is the same scheme the eight-fault pipeline uses; the
  smoke test verifies it works for one.
- **Solver-side wiring is out of scope.** The smoke test ends at "the existing BP5
  mesh loader reads the file without error" — manual confirmation only.

## Phase 1: Wire up `--fault` filter on the pipeline scripts

### Goal
After this phase the four pipeline driver scripts (`audit_ts_quality.py`,
`ts_to_stl.py`, `generate_safs_mesh.py`, `validate_msh.py`, plus
`write_fault_provenance.py` from `PLAN_domain.md` Phase 3) accept a `--fault NAME`
or `--include-fault NAME [NAME ...]` flag that restricts the run to a single CFM
fault by short-name (e.g., `safs_sbmt_millcreek`).

### Files to Modify
- `mesh/audit_ts_quality.py` — accept `--include-fault` (repeatable). When given, the
  audit walks only the matching `.ts` files; the output CSV still has the same schema.
- `mesh/ts_to_stl.py` — same flag. Skips writing STLs for excluded faults; updates
  `mesh/transform.json` `extra.input_files` to reflect the actual subset; writes
  `mesh/bbox.json` only for the included faults.
- `mesh/generate_safs_mesh.py` — same flag, plus `--out-suffix` to control the
  output mesh filename and the output directory (`output/smoke_millcreek/...`).
- `mesh/safs.geo` — no source change. The `.geo` reads which STLs to merge from a
  comma-separated `-setstring` parameter that the driver injects, so the driver
  controls which faults are geometric inputs.
- `mesh/write_fault_provenance.py` (defined in `PLAN_domain.md`) — gains the same
  `--include-fault` flag for symmetry.
- `mesh/validate_msh.py` — gains an `--expected-faults NAME [...]` flag listing which
  short-names should appear in `fault_provenance.json`. Default = all eight; smoke
  test passes `--expected-faults safs_sbmt_millcreek`.

### Detailed Requirements

1. **Short-name table is the canonical fault identifier.** Defined once in
   `mesh/safs_origin.py` (per `PLAN_origin.md` Phase 1), reused everywhere:
   ```python
   FAULT_SHORT_NAMES = {
       "SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6":              "safs_mjvs_saf",
       "ETRA-PMFZ-MULT-Pinto_Mountain_fault-CFM5":           "safs_pmfz_pinto",
       "SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4":        "safs_sbmt_millcreek",
       "SAFS-SAFZ-SBMT-Mission_Creek_fault_strand-CFM4":     "safs_sbmt_missioncreek",
       "SAFS-SAFZ-COAV-Mission_Creek_fault_strand-CFM4":     "safs_coav_missioncreek",
       "SAFS-SAFZ-MULT-Banning_fault-CFM6":                  "safs_mult_banning",
       "SAFS-SAFZ-MULT-Southern_San_Andreas_fault_and_Banning-CFM6":
                                                              "safs_mult_ssaf_banning",
       "SAFS-SAFZ-SBMT-San_Andreas_fault-CFM6":              "safs_sbmt_saf",
   }
   ```
   `--include-fault` accepts the short-name (NOT the CFM ID) for ergonomics.

2. **`safs.geo` STL list parameter.** The `.geo` accepts a single string (Gmsh's
   `DefineConstant` with `string` type, or driver-side templating — Gmsh `-setstring`
   on top-level strings is supported in 4.10+):
   ```
   DefineConstant[ stl_list = {"safs_mjvs_saf,safs_pmfz_pinto,safs_sbmt_millcreek,...",
                               Name "Comma-separated STL short-names"} ];
   ```
   Inside `safs.geo`, split the string into tokens and `Merge` each one. (Gmsh's geo
   language has limited string ops; the driver will instead emit a small `safs_includes.geo`
   side-file containing one `Merge "stl/<short>.stl";` per included fault, and the main
   `safs.geo` does `Include "safs_includes.geo";`. Cleaner than parsing strings inside `.geo`.)

3. **Driver flow with `--include-fault`.**
   ```python
   # generate_safs_mesh.py, sketch
   included = parse_include_fault(args.include_fault)            # list[str], short names
   bbox     = read_bbox_json("mesh/bbox.json")                   # dict[short→bbox]
   bbox     = {k: v for k, v in bbox.items() if k in included}   # filter
   # compute domain box from the FILTERED bbox (per PLAN_domain.md Phase 1)
   write_includes_geo("mesh/safs_includes.geo", included)        # one Merge per short name
   # then invoke gmsh as before
   ```
   The same `--include-fault` set must have been used in the prior `ts_to_stl.py` run,
   otherwise some `stl/<short>.stl` files won't exist. Hard-fail the driver before
   invoking Gmsh:
   ```python
   for s in included:
       p = stl_dir / f"{s}.stl"
       if not p.exists():
           sys.exit(f"missing STL {p}; rerun ts_to_stl.py --include-fault {s}")
   ```

4. **Backwards compatibility.** When `--include-fault` is omitted, all eight faults
   are processed (default unchanged from the master plan).

### Acceptance Criteria
- [ ] `python audit_ts_quality.py --include-fault safs_sbmt_millcreek --res 2000`
      writes a CSV with exactly one row.
- [ ] `python ts_to_stl.py --include-fault safs_sbmt_millcreek --res 2000` writes
      exactly `mesh/stl/safs_sbmt_millcreek.stl` and updates `mesh/bbox.json` and
      `mesh/transform.json` accordingly.
- [ ] `python generate_safs_mesh.py --include-fault safs_sbmt_millcreek
      --out-suffix smoke_millcreek` produces `output/smoke_millcreek/safs_2000m.msh`.
- [ ] Without any `--include-fault` flag, every script behaves exactly as in the
      master plan (regression check via the existing eight-fault acceptance criteria).

### Dependencies
- Depends on: `PLAN.md` Phases 1, 2; `PLAN_origin.md` Phases 1–3;
  `PLAN_domain.md` Phases 1–3.
- Required by: Phase 2 of this document.

---

## Phase 2: Run the smoke test and capture artefacts

### Goal
After this phase the directory `mesh/output/smoke_millcreek/` contains every artefact
the multi-fault pipeline would produce, but for the single Mill Creek strand. The run
is reproducible from a single shell command.

### Files to Create
- `mesh/run_smoke_millcreek.sh` — one-shot wrapper that invokes the four scripts in
  order with the right flags. Idempotent (safe to re-run).

### Detailed Requirements

1. **Wrapper contents.**
   ```bash
   #!/usr/bin/env bash
   set -euo pipefail
   cd "$(dirname "$0")"
   FAULT=safs_sbmt_millcreek
   RES=2000
   OUTDIR=output/smoke_millcreek
   mkdir -p "$OUTDIR"

   python audit_ts_quality.py \
       --cfm-dir "$HOME/Documents/Earthquake Cycle Modeling of San Andreas Fault System/CFM_data" \
       --res "$RES" \
       --include-fault "$FAULT" \
       --out "$OUTDIR/cfm_audit_${RES}m.csv"

   python ts_to_stl.py \
       --cfm-dir "$HOME/Documents/Earthquake Cycle Modeling of San Andreas Fault System/CFM_data" \
       --res "$RES" \
       --include-fault "$FAULT" \
       --out-dir "$OUTDIR/stl"

   python generate_safs_mesh.py \
       --include-fault "$FAULT" \
       --stl-dir "$OUTDIR/stl" \
       --bbox-json "$OUTDIR/bbox.json" \
       --transform-json "$OUTDIR/transform.json" \
       --res-f 1000 \
       --res-ff 20000 \
       --ramp-dist 30000 \
       --buf-x 50000 --buf-y 50000 --depth 50000 \
       -o "$OUTDIR/safs_smoke_2000m.msh"

   python write_fault_provenance.py \
       --msh "$OUTDIR/safs_smoke_2000m.msh" \
       --stl-dir "$OUTDIR/stl" \
       --include-fault "$FAULT" \
       --out "$OUTDIR/fault_provenance.json"

   python validate_msh.py \
       --msh "$OUTDIR/safs_smoke_2000m.msh" \
       --transform-json "$OUTDIR/transform.json" \
       --provenance-json "$OUTDIR/fault_provenance.json" \
       --expected-faults "$FAULT" \
       --report "$OUTDIR/validation_report.txt"
   ```

2. **Environment.** Per `miniapps/seas/CLAUDE.md`: `conda activate pythonenv` for Gmsh
   and the Python tooling. The wrapper does not activate the env (the user does);
   it will fail loudly if `gmsh` is not on PATH.

3. **Reproducibility seal.** `validate_msh.py` writes
   `output/smoke_millcreek/safs_smoke_2000m.provenance.json` per `PLAN.md` Phase 6,
   recording git commit, Gmsh version, input SHA-256, and output SHA-256. The smoke
   test passes if a second run produces a `.msh` with the same vertex count and
   tag inventory (HXT non-determinism allows ±0.5% tet count drift).

### Acceptance Criteria
- [ ] `bash mesh/run_smoke_millcreek.sh` exits 0 on a clean checkout.
- [ ] The five artefacts exist in `output/smoke_millcreek/`:
      `cfm_audit_2000m.csv`, `stl/safs_sbmt_millcreek.stl`, `safs_smoke_2000m.msh`,
      `fault_provenance.json`, `validation_report.txt`.
- [ ] Re-running the wrapper without deleting outputs is a no-op for the first script
      (idempotent audit) and a deterministic re-run for the rest (overwrites).

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 3.

---

## Phase 3: Targeted invariant checks specific to a single fault

### Goal
After this phase, `validate_msh.py --expected-faults safs_sbmt_millcreek` runs a
suite of checks tuned for the single-fault scenario and emits a green/red report. The
checks below are *in addition to* the generic acceptance criteria already in
`PLAN_domain.md` Phase 3.

### Files to Modify
- `mesh/validate_msh.py` — add a `--single-fault-mode` block that runs the checks
  in (1)–(8) below when `--expected-faults` has length 1.

### Detailed Requirements — checks the smoke test must run

1. **Tag inventory.**
   - `.msh` physical groups exactly `{1, 2, 3, 4, 5, 6, 10, 100}`.
   - No tag is empty (each lateral face has > 0 triangles, ztop and zbot > 0
     triangles, fault has > 0 triangles, volume has > 0 tets).

2. **Domain box arithmetic.**
   - Read `domain_box.json`. With Mill Creek's local bbox
     `X∈[-30 212, +28 501], Y∈[+2 227, +19 534]` and the default 50 km buffer / 10 km
     snap, the smoke-test box should be:
     ```
     X ∈ [-90 000, +90 000] m
     Y ∈ [-50 000, +70 000] m
     Z ∈ [-50 000,       0] m
     ```
     (`floor((-30 212 - 50 000)/10 000)*10 000 = -90 000`, etc.)
   - Verify `domain_box.json` matches this exactly.

3. **Fault clearance.**
   - Mill Creek's bounding box must sit ≥ 5 km clear of every lateral face and z_bot.
     With the box above, X clearance ≥ 60 km, Y clearance from `ymin = -50 km` to
     fault `ymin = +2 227 m` is 52 km, Z clearance from z_bot to fault zmin is 32 km.
     All comfortably clear. Hard-fail otherwise.

4. **Free-surface trace exists.**
   - Mill Creek dips at 83°, so its surface trace is ~at z = 0. Verify that
     **at least one tag-100 triangle has all three vertices with z > -100 m**
     (i.e., the fault reaches the free surface). This catches a failure mode where
     the z>0 clamp accidentally pushes the trace below z=0 by some tolerance.
   - Verify that **no tag-100 triangle has any vertex with z > 0** (clamp worked).

5. **Fault is internal (the SEAS critical invariant).**
   - For every tag-100 triangle T: T must have exactly 2 adjacent tetrahedra
     (one on each side). A triangle with only 1 adjacent tet means the mesher made
     it a boundary, which would be solver-fatal. Implementation: build face→tet
     adjacency from the `.msh` element list and assert every fault face has 2 tets.
   - Exception: triangles whose all three vertices have z = 0 (the trace) lie on
     `ztop` and may legitimately have only 1 tet. Count these separately and
     report. Mill Creek's trace will produce O(50) such trace triangles given its
     ~60 km along-strike length and 1 km mesh size.

6. **Mesh sizing on the fault.**
   - For tag-100 triangles, mean edge length should be `1000 ± 300 m` (per
     `PLAN_domain.md` Phase 2 acceptance). Histogram-print the distribution.
   - Max edge length on the fault should be `< 2 000 m`. A larger value indicates
     the size field did not bind on the fault — likely a bug in `Field[1].SurfacesList`.

7. **Mesh sizing far from fault.**
   - For tetrahedra with all 4 vertices > 30 km from any tag-100 face, mean edge
     length should be `> 12 000 m`. Catches the case where the threshold field
     never reached `res_ff`.

8. **Provenance partition.**
   - `fault_provenance.json["faults"]["safs_sbmt_millcreek"]["n_triangles_in_msh"]`
     == total number of tag-100 triangles. (Trivial when only one fault is included,
     but still a useful regression for the partition logic.)
   - The list of `triangle_indices_in_msh` is sorted, no duplicates, and indexes
     valid triangles in the `.msh`.

9. **Round-trip UTM recovery.**
   - Pick 5 random tag-100 triangle vertices. Apply `local_to_utm` from `safs_origin`
     using `transform.json`. Confirm the recovered UTM coordinate lies inside the
     original CFM bbox (`X ∈ [469 788, 528 501], Y ∈ [3 767 227, 3 784 534]`,
     `Z ∈ [-18 014, 2 593]`) with no margin. Catches an origin-application sign bug.

10. **Tet quality.**
    - `min Gamma >= 0.10`, `mean Gamma >= 0.55` (per `PLAN_domain.md`). Mill Creek's
      near-vertical dip and ~60 km × 18 km extent are the *easiest* single-fault
      stress test for HXT — failure here points at infrastructure rather than CFM
      complexity.

### Implementation
`validate_msh.py` accumulates each check's pass/fail into a list, prints a summary
table, and exits 0 only if all pass. Failed checks include the offending values
(e.g., "check 5 FAIL: 12 fault triangles have only 1 adjacent tet, indices [...]").

### Acceptance Criteria
- [ ] All ten checks above pass on the smoke-test mesh.
- [ ] `validation_report.txt` contains the full check-by-check table with a green
      summary line at the end.
- [ ] Deliberately breaking each check (e.g., setting `--res-f 50000` to violate
      check 6) flips that check to FAIL with a useful message — verifies the
      checks are not vacuously passing.

### Dependencies
- Depends on: Phases 1, 2.
- Required by: nothing.

---

## Phase 4: Catalogue of expected bug classes (what the smoke test is hunting)

This phase is documentation only — no code, no acceptance criteria. The list of
checks in Phase 3 is designed to surface specific, named failure modes. For each, we
record: **(a) what would cause it, (b) which check catches it, (c) where to fix it.**

| # | Bug class                                       | Caught by check  | Fix lives in                                                     |
| - | ----------------------------------------------- | ---------------- | ---------------------------------------------------------------- |
| 1 | Origin sign error (subtract vs add)             | 9                | `safs_origin.py:utm_to_local`                                    |
| 2 | Z translation accidentally applied              | 4                | `safs_origin.py` — `Z0 = 0` literal                              |
| 3 | z > 0 clamp not applied (or applied to wrong axis) | 4              | `ts_to_stl.py` — clamp ordering                                  |
| 4 | STL vertex deduplication breaks triangles       | Phase 2 STL audit | `ts_to_stl.py` — degenerate-after-collapse handling             |
| 5 | `Distance.SurfacesList` empty (fragmentation collapsed all tags) | 6 |  `safs.geo` — surface tag recovery after fragments         |
| 6 | `Threshold.SizeMin/SizeMax` swapped             | 6 / 7            | `safs.geo` — Threshold field constants                           |
| 7 | `Background Field` not set                      | 6                | `safs.geo` — `Background Field = 2;`                             |
| 8 | Box face bbox query catches a fault triangle    | 1 / 5            | `safs.geo` — `eps`/`round_to` discipline                         |
| 9 | Fault becomes a boundary (1 adjacent tet)       | 5                | `safs.geo` — `BooleanFragments` order or fault embed             |
| 10 | Fault tag != 100                               | 1                | `safs.geo` — `Physical Surface("fault", 100)` literal            |
| 11 | Multiple physical tags per element confusion   | 1                | `safs.geo` — only one `Physical Surface` line per surface set    |
| 12 | `transform.json` not written / wrong schema    | 9                | `ts_to_stl.py:write_transform_json`                              |
| 13 | Provenance partition leaks (sum != total)      | 8                | `write_fault_provenance.py:_assign_to_fault`                     |
| 14 | Driver/`.geo` mismatch on box extents          | 2                | `generate_safs_mesh.py:_compute_box`                             |
| 15 | OCC tolerance too coarse (fault glued to box)  | 5                | `safs.geo` — `Geometry.Tolerance = 1.0`                          |
| 16 | Free-surface trace lost (clamped below 0)      | 4                | `ts_to_stl.py` — clamp threshold                                 |
| 17 | Mesh size unit confusion (km vs m)             | 6                | `safs.geo` / driver — units docstring                            |

The smoke test is considered successful only when every bug in this catalogue
either does not occur in the Mill Creek run, or is detected by the listed check.

---

## Testing Strategy
- **Pre-flight:** before running the wrapper, manually verify the input file exists at
  `~/Documents/.../CFM_data/SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4_2000m.ts` and
  that `gmsh -version` reports ≥ 4.10.
- **Run:** `bash mesh/run_smoke_millcreek.sh`. Total wall time should be < 2 minutes
  for the dev tier (1 km on-fault, 20 km far-field). Anything above 5 minutes
  indicates a sizing-field mistake — abort and inspect.
- **Reproducibility check:** run twice; compare `provenance.json` vertex/element
  counts; allow ±0.5% drift due to HXT non-determinism.
- **Adversarial sanity:** intentionally pass `--res-f 50000` (oversized) and confirm
  that check 6 fails. Restore default after.

## Risk Assessment
1. **Mill Creek's z>0 vertices are 10% of the surface** — the largest fraction of all
   eight CFM faults. If the z>0 clamp has a bug, this fault is the most likely to
   expose it (versus, say, MJVS-SAF where z>0 is ≪ 1%). Hence Mill Creek as the
   smoke fault, not just because of its medium triangle count.
2. **Mill Creek dips at 83° (near-vertical).** Its triangulation is well-conditioned
   — the surface lies almost in a plane. If something fails for Mill Creek, the
   problem is in the *infrastructure*, not in handling steeply-curved freeform
   surfaces. This makes the smoke test a clean infrastructure check before we
   complicate it with the SBMT-SAF (51° dip) and the San Gorgonio knot.
3. **The `--include-fault` flag plumbing touches every script.** The most likely bug
   surface in this plan is *not* the meshing itself but the per-script flag handling
   — a flag that filters STL writes but not bbox writes will desynchronize the
   pipeline. The Phase 1 acceptance criteria specifically test this.
4. **A passing smoke test does NOT prove the multi-fault pipeline works.** Mill
   Creek alone never exercises `BooleanFragments` between two CFM surfaces, never
   exercises shared edges between strands, never exercises the eastern half of the
   domain. After the smoke test passes, the next test should be a two-fault run
   (Mill Creek + SBMT-SAF) to exercise pairwise intersection. That two-fault run
   is **out of scope** of this document.
