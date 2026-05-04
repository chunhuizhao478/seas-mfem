# SAFS Realistic-Geometry Mesh Pipeline

3-D conformal tetrahedral finite-element meshes of the San Andreas Fault
System, built from SCEC Community Fault Model (CFM) GOCAD-TSurf surfaces
using Gmsh.  Inputs: 8 CFM `.ts` files at 500m / 1000m / 2000m
resolutions.  Outputs: Gmsh-2.2 `.msh` files with a BP5-compatible tag
scheme, plus ParaView `.vtu` and Dolfin `.xml` for visualization.

This README covers (1) what's working today, (2) how to run it, (3) the
limits we hit, and (4) the roadmap for the missing piece (multi-fault
conformal mesh).

---

## TL;DR

```bash
conda activate pythonenv
cd miniapps/seas/safs/mesh

# Single fault — passes all 10 SEAS-grade validation checks.
bash run_smoke_millcreek.sh
# → output/smoke_millcreek/output/safs_smoke_1000m.{msh,vtu,xml}

# All 8 CFM faults in one mesh — viewable, but fault triangles are NOT
# proper internal interfaces (see "Caveats" below).
bash run_all8_2000m.sh
# → output/all8_2000m/output/safs_all8_2000m.{msh,vtu,xml}

# Spatially-disjoint subsets — each subset is SEAS-grade, but only 6 of
# 8 CFM faults are covered (Pinto + Banning excluded as CFM-defective).
bash run_subsets_2000m.sh
# → output/subset_{NW,South,Mill,MissionSBMT}_2000m/output/...
```

To open in ParaView and colour the fault by tag:
`paraview <path>.vtu`, then in the Properties panel set Coloring to
`physical`.  Tags: 1–6 = box faces, 10 = bulk volume, 100 = fault.

---

## What works today, by mode

| wrapper | scope | SEAS-ready? | gmsh 3D algo | notes |
|---|---|---|---|---|
| `run_smoke_millcreek.sh`  | 1 fault (Mill Creek strand, 1000m) | **YES** (10/10 validator checks) | HXT | the trusted reference |
| `run_subsets_2000m.sh`    | 6 faults across 4 disjoint subsets, 2000m | **YES per subset** | HXT | drops Pinto + Banning entirely (CFM source defects); 4 separate `.msh` files |
| `run_all8_2000m.sh`       | all 8 faults in one mesh, 2000m | **NO** — fault tris are not internal interfaces | Frontal-Delaunay | for visualization / scoping; needs splitting work to upgrade to SEAS-ready |

"SEAS-ready" means: every fault triangle has exactly two adjacent tets
(one on each side), so a DG friction solver can compute traction jumps
across the fault.  This is validator check 5 (internal-interface
invariant).

---

## Pipeline architecture

Each wrapper invokes the same six-stage pipeline, with different
`--include-fault` flags and resolution:

```
   audit_ts_quality.py        per-CFM-file QC + overlap-pair detector
            │
            ▼
   ts_to_stl.py               clean → utm_to_local → clamp z → dedup →
                              drop degenerate → repair coplanar overlaps →
                              write per-fault ASCII STL + transform.json
            │
            ▼
   generate_safs_mesh.py      compute domain box (snap to 10 km grid,
                              50 km buffers); when len(faults) > 1, pre-
                              merge per-fault STLs into safs_combined.stl
                              with shared-vertex welding at 1 cm; emit
                              safs_includes.geo; invoke gmsh -3 safs.geo
            │
            ▼
   safs.geo                   Box(1) bulk volume + Embed fault surfaces
                              + Distance/Threshold size-field
                              (res_f near fault, res_ff far field) +
                              BP5 physical groups
            │
            ▼
   write_fault_provenance.py  KDTree-assign every tag-100 triangle to its
                              source CFM fault; output schema-versioned
                              fault_provenance.json
            │
            ▼
   validate_msh.py            10-check suite (tag inventory, box arithmetic,
                              clearance, free-surface trace, interface
                              invariant, edge sizing, far-field edge,
                              provenance partition, UTM round-trip,
                              tet quality)
            │
            ▼
   convert_msh.py             Gmsh .msh  →  ParaView .vtu  + Dolfin .xml
```

Single source of truth for the local frame: `safs_origin.py` defines
`E0=500_000, N0=3_765_000, Z0=0` (UTM Zone 11N false easting; near
geographic centroid of the 8-fault complex).  See `PLAN_origin.md`.

Single source of truth for the BP5 tag scheme: `safs.geo` Physical-block
matches `bp5/mesh/bp5_v2.geo` exactly (1–6 boundary, 10 volume, 100 fault).

Single source of truth for fault short-names: `FAULT_SHORT_NAMES` dict in
`safs_origin.py` — 8 CFM IDs ↔ 8 short-names used by every script.

---

## How to run

### Prerequisites

```bash
conda activate pythonenv     # provides gmsh ≥ 4.10, meshio, scipy,
                              # jsonschema, numpy
```

CFM data must be at
`~/Documents/Earthquake Cycle Modeling of San Andreas Fault System/CFM_data/`
(override via `SAFS_CFM_DIR=<path>` env var in any wrapper).

### Smoke test (recommended starting point)

```bash
bash run_smoke_millcreek.sh
# Expected: "10/10 checks passed" near the end.
# Run time: ~10 s on a laptop.
```

### All 8 faults in one mesh (visualization)

```bash
bash run_all8_2000m.sh
# Expected: "8/10 checks passed".  Checks 5 and 7 fail (see Caveats).
# Run time: ~5 s.
```

### Per-subset SEAS-ready meshes (six of eight faults)

```bash
bash run_subsets_2000m.sh
# Expected: "10/10 checks passed" for each of the 4 subsets.
# Run time: ~30 s total.
```

### Custom subset

```bash
python audit_ts_quality.py    --cfm-dir "$CFM_DIR" --res 1000 \
    --include-fault safs_mjvs_saf --include-fault safs_sbmt_saf \
    --out my_audit.csv
python ts_to_stl.py           --cfm-dir "$CFM_DIR" --res 1000 \
    --include-fault safs_mjvs_saf --include-fault safs_sbmt_saf \
    --out-dir my_stl --transform-json my_transform.json \
    --bbox-json my_bbox.json --free-surface-clearance 100 \
    --repair-overlaps
python generate_safs_mesh.py  --include-fault safs_mjvs_saf \
    --include-fault safs_sbmt_saf --stl-dir my_stl \
    --bbox-json my_bbox.json --transform-json my_transform.json \
    --res-f 1000 --res-ff 20000 --ramp-dist 30000 \
    -o my_output.msh
python write_fault_provenance.py --msh my_output.msh \
    --stl-dir my_stl --transform-json my_transform.json \
    --include-fault safs_mjvs_saf --include-fault safs_sbmt_saf \
    --out my_provenance.json
python validate_msh.py        --msh my_output.msh \
    --transform-json my_transform.json --bbox-json my_bbox.json \
    --provenance-json my_provenance.json
python convert_msh.py         --msh my_output.msh
```

Every script supports `--help`.

### Tunables (env var → wrapper default)

| env var          | default | meaning |
|---|---|---|
| `SAFS_RES`        | 1000 / 2000 | CFM resolution suffix in metres (500m, 1000m, 2000m available) |
| `SAFS_CLEARANCE`  | 100     | fault-trace stand-off below z=0 (m); HXT needs > 0 |
| `SAFS_RES_F`      | 1000 / 2000 | target on-fault edge length (m) |
| `SAFS_RES_FF`     | 20000 / 25000 | far-field edge length (m) |
| `SAFS_RAMP`       | 30000 / 35000 | size-field linear ramp width (m) |
| `SAFS_BUF`        | 50000   | lateral buffer around fault footprint (m) |
| `SAFS_DEPTH`      | 50000   | bulk box depth (m) |
| `SAFS_CFM_DIR`    | `~/.../CFM_data` | CFM source directory |

### Eight CFM faults (the input set)

```
SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6                 → safs_mjvs_saf
ETRA-PMFZ-MULT-Pinto_Mountain_fault-CFM5              → safs_pmfz_pinto
SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4           → safs_sbmt_millcreek
SAFS-SAFZ-SBMT-Mission_Creek_fault_strand-CFM4        → safs_sbmt_missioncreek
SAFS-SAFZ-COAV-Mission_Creek_fault_strand-CFM4        → safs_coav_missioncreek
SAFS-SAFZ-MULT-Banning_fault-CFM6                     → safs_mult_banning
SAFS-SAFZ-MULT-Southern_San_Andreas_fault_and_Banning-CFM6 → safs_mult_ssaf_banning
SAFS-SAFZ-SBMT-San_Andreas_fault-CFM6                 → safs_sbmt_saf
```

---

## Caveats — what's NOT working today

### 1. Multi-fault internal-interface conformity

`run_all8_2000m.sh` produces a complete mesh, but Frontal-Delaunay
(Algorithm3D=4) does NOT make the fault triangles shared faces of the
bulk tets.  Validator check 5 fails: every tag-100 triangle has 0 (not
2) adjacent tets.  Visualization works; the SEAS DG solver cannot
compute traction jumps across the fault as-is.

The HXT-based path (Algorithm3D=10), which DOES enforce conformity, was
abandoned for multi-fault input because HXT rejects every multi-fault
combination with PLC errors:

- 10 fault-pair combinations have 3-D triangle-triangle crossings:
  Mill Creek × SBMT-SAF (163 pairs), Mission Creek SBMT × SAF (135),
  Pinto × Mill Creek (42), SBMT-Mill × SBMT-MissionCreek (41), etc.
- 4 fault-pair combinations have transverse trace crossings at the
  free-surface clearance plane.
- A simplified "drop the lower-rank intersector" trim was implemented
  and discarded — even at 17.3% triangle drop it could not eliminate
  the T-junctions the trim itself produces on the cut boundaries.

### 2. Two CFM faults are CFM-source-defective

`safs_pmfz_pinto` (Pinto Mountain) and `safs_mult_banning` (Banning)
fail HXT PLC recovery on their **own** at every CFM resolution
(500m, 1000m, 2000m).  The defect is in the source `.ts` triangulation:
multiple coplanar same-side overlapping triangle pairs and intra-fault
trace self-crossings beyond what `--repair-overlaps` removes.  The
SSAF+Banning composite surface is fine and is included in the
all-8 build.

### 3. Tag-100 fault tris in `run_all8_2000m.sh` may have slivers

Frontal-Delaunay produces a slightly worse mesh in the fault-fault
crossing zones than HXT does for clean inputs, but check 10 confirms
0 sliver tets in the all-8 build (mean γ=0.87, min γ=0.54) — ironically
better than the smoke test.  Tet quality is fine; what's missing is
*topological* fault-tet adjacency.

---

## Roadmap

The piece missing today is a **proper multi-fault conformal mesh** —
all 8 CFM faults in one mesh, with each fault triangle a true internal
interface (validator check 5 passes).  Three paths, ranked by effort
and risk.

### Path A — Triangle-splitting preprocessor ("Option 3 done right")

For each crossing fault pair (L, H):
1. compute the 3-D intersection segment (where the planes of the two
   triangles cross),
2. split each triangle in BOTH faults at that segment, producing
   2-3 conformal sub-triangles per side,
3. weld the new edges into a shared-vertex curve so the two faults'
   triangulations agree exactly along the seam,
4. cascade: a triangle split by intersection with H may need further
   splitting against M, W, etc.,
5. repair degenerate sub-triangles produced by near-coplanar splits.

The only geometric primitive that survives in this repo is
`audit_ts_quality.py:find_overlap_pairs`.  An earlier `cfm_trim.py`
contained `tri_tri_intersect_3d`, `edge_pierces_tri`,
`_point_to_tri_distance`, and `seg_seg_cross_2d`, but the entire
`miniapps/seas/safs/` tree is currently untracked in git, so those
helpers are NOT recoverable and would have to be rewritten.

**Estimated effort:** 1–2 weeks for a robust first version with unit
tests on synthetic two-plane configurations and integration on real
CFM data.  **Numerical-robustness risk:** triangle-triangle
intersection is fragile in floating point when triangles are nearly
coplanar (which is exactly the SAFS case — adjacent CFM pieces are
both ~80–90° dip).  Exact arithmetic via `gmpy2.mpfr` predicates may
be required.

**Required validation:** all 10 validator checks pass on the all-8
mesh; provenance partition lossless across splits.

### Path B — External boolean preprocessor

Use `manifold` (Google's new library), `PyMesh`/CGAL, or `Cork` to
compute the surface union of the 8 STLs with conformal intersection
curves built in.  This replaces ~500 lines of custom carving code with
a 5-line `import manifold; combined = sum(faults)`.

**Estimated effort:** 2–3 days, dominated by macOS install pain
(`manifold` needs a recent CMake + a C++17 compiler; PyMesh has
multiple system-library dependencies).  **Risk:** these tools assume
**closed manifolds** for boolean union; CFM fault surfaces are open
(have boundary edges).  Some preprocessing may be needed to virtually
"close" each surface before the union, then reopen the seam.  Worth
investigating but not a guaranteed save.

**Required validation:** same as Path A.

### Path C — Different mesher

TetGen with PLC has different topology requirements than Gmsh-HXT and
may tolerate some configurations Gmsh rejects.  Or Coreform's
commercial SCEC-CFM workflow.  Or a fault-mortar approach where each
fault is meshed separately and the SEAS solver couples them via a
mortar method (changes the solver as well as the mesher).

**Estimated effort:** TetGen swap-in: 1–2 days to wire `tetgen` calls
through `generate_safs_mesh.py` and re-run validation.  Mortar is a
large solver-side project, out of scope here.

**Risk:** TetGen's PLC documentation suggests it has the same
non-self-intersection requirement as HXT.  Limited upside.

### Recommendation

Start with **Path B** (`manifold`) — lowest code investment, leverages
mature tooling.  If the open-surface boundary problem can't be
finessed, fall back to **Path A** (custom triangle splitter) with the
existing primitives.

---

## File map

```
miniapps/seas/safs/mesh/
├── README.md                      ← you are here
├── PLAN.md, PLAN_origin.md,       ← design documents in ../
│   PLAN_domain.md,
│   PLAN_smoke_millcreek.md
│
├── safs_origin.py                 ← origin module (R-001..R-009 baseline)
├── transform.schema.json          ← JSON Schema for transform.json
├── audit_ts_quality.py            ← TSurf parser + QC + overlap detector
├── ts_to_stl.py                   ← cleanup → STL with --repair-overlaps
├── safs.geo                       ← Gmsh discrete-Embed + Algorithm3D=4
├── generate_safs_mesh.py          ← driver: box, sizing, multi-STL combine
├── write_fault_provenance.py      ← KDTree-based per-tri attribution
├── validate_msh.py                ← 10-check post-mesh validator
├── convert_msh.py                 ← .msh → .vtu + .xml
│
├── run_smoke_millcreek.sh         ← 1-fault SEAS-grade smoke test
├── run_all8_2000m.sh              ← 8-fault visualization mesh
├── run_subsets_2000m.sh           ← 4 SEAS-grade subset meshes (6/8 faults)
│
└── output/                        ← (generated) per-run outputs
    ├── smoke_millcreek/
    ├── all8_2000m/
    └── subset_{NW,South,Mill,MissionSBMT}_2000m/
```

---

## References

- SCEC CFM data portal: https://www.scec.org/research/cfm
- BP5 SCEC SEAS spec: https://strike.scec.org/cvws/seas/download/SEAS_BP5_QD.pdf
- Project SEAS-MFEM CLAUDE.md: `miniapps/seas/CLAUDE.md`
- Tandem reference: `/Users/chunhuizhao/projects/tandem/`
- Code review and roadmap: `REVIEW.md` at repo root
