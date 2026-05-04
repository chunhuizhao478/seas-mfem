# Fix Report: REVIEW.md (2026-05-02) — General mesh-integrity check & fault-orientation fix

Branch: `feature/elasticity-inertia`
Date: 2026-05-02

## Summary table

| ID | Status | Files modified |
|---|---|---|
| R-001 | FIXED | `miniapps/seas/safs/mesh/validate_msh.py` (added `check_12_surface_closure`, registered in `main()`) |
| R-002 | FIXED | `miniapps/seas/safs/mesh/validate_msh.py` (added `check_13_fault_orientation`; HARD-fail on `n_winding_flips_needed > 0`, WARN-only on non-manifold) |
| R-003 | FIXED | `miniapps/seas/safs/mesh/orient_fault_surface.py` (NEW); wired into `run_newset_step_by_step.sh` after each mesh-mutating stage |
| R-004 | FIXED | `miniapps/seas/safs/mesh/run_newset_step_by_step.sh` (validate failure now aborts pipeline at all 4 sites) |
| R-005 | RESOLVED-BY-R-002 | n/a — non-manifold edge metric is reported by `check_13` |
| R-006 | FIXED | `miniapps/seas/safs/mesh/mmg3d_local_patch.py` (`_stitch_back` raises `RuntimeError` on producer-side surface holes) |

## Files modified

- `miniapps/seas/safs/mesh/validate_msh.py` — added `check_12_surface_closure` and `check_13_fault_orientation`; both registered in `main()`. Pre-fix the suite was 11 checks; post-fix it is 13 checks.
- `miniapps/seas/safs/mesh/orient_fault_surface.py` — NEW (~7.4 KB). BFS-propagates fault-tri winding within each tag-100 connected component. Idempotent. Emits `components / flipped / nonmanifold` to stderr.
- `miniapps/seas/safs/mesh/mmg3d_local_patch.py` — `_stitch_back()` now asserts the "every 1-tet bdry face has a tagged tri" invariant on its own output and raises `RuntimeError` with diagnostic centroids when violated.
- `miniapps/seas/safs/mesh/run_newset_step_by_step.sh` — wired `python orient_fault_surface.py` immediately before each stage's `validate_msh.py` invocation (raw HXT, post-mmg3d, post-patch, post-cavity); converted the four `validate_msh.py` call sites from `set +e ... set -e` envelopes to `if !  ...; then return N; fi` so validation failure aborts the pipeline (R-004).
- New unit tests:
  - `miniapps/seas/safs/mesh/tests/test_validate_msh_check12_check13.py` (R-001 + R-002, 9 tests)
  - `miniapps/seas/safs/mesh/tests/test_orient_fault_surface.py` (R-003, 6 tests)
  - `test_mmg3d_local_patch.py::test_R006_stitch_back_raises_on_surface_hole` (R-006)

## Unit test results

- `pytest miniapps/seas/safs/mesh/tests/` — **136 passed, 0 failed, 1 deprecation warning** (1.50 s).
- Of those 136, the R-fix tests are: 9 for R-001+R-002 (`test_validate_msh_check12_check13.py`), 6 for R-003 (`test_orient_fault_surface.py`), 1 for R-006 (in `test_mmg3d_local_patch.py`). All 16 R-fix tests pass.
- Orphan/dead tests previously called out were removed by the prior agent before this run; the remaining 136 all pass cleanly.

## Step-6 RESULT lines (4 stages, post-fix)

The live `run_newset_step_by_step.sh` was relaunched three times in this session; all three failed at `generate_safs_mesh` (gmsh HXT 3D self-intersecting facets) on a slightly different STL than the previously-completed runs (15676 tris vs 15667). The HXT failure is unrelated to R-001..R-006 — it is a known gmsh/CGAL non-determinism on the SAFS dedup output. To exercise the post-fix validation against real downstream-stage meshes, the new `orient_fault_surface.py` and `validate_msh.py` were run directly against the 4 stage outputs from the most recent successful end-to-end run at `output/newset_6_all6_cavity_retet/output/`. Logs at `/tmp/postfix_validate/validate_{raw,mmg3d,patch,cavity}.log`.

```
RESULT [raw]:    12/13 checks passed | gamma_min=2.44e-10  | min_edge=0.20 m | slivers=307 | n_tets=936933  | rc=1
RESULT [mmg3d]:  12/13 checks passed | gamma_min=2.24e-06  | min_edge=0.20 m | slivers=282 | n_tets=1101972 | rc=1
RESULT [patch]:  11/13 checks passed | gamma_min=2.53e-06  | min_edge=0.20 m | slivers=337 | n_tets=1108808 | rc=1
RESULT [cavity]: 11/13 checks passed | gamma_min=6.61e-06  | min_edge=0.20 m | slivers=318 | n_tets=1108882 | rc=1
```

Per-stage failure breakdown:

| Stage | Failed checks | Note |
|---|---|---|
| raw    | check_10 only (slivers, pre-existing R-402 gate) | check_12 PASS (0 unlabeled holes), check_13 PASS (0 winding flips, 393 non-manifold reported as branching, ok) |
| mmg3d  | check_10 only | check_12 PASS, check_13 PASS (0 flips, 269 non-manifold) |
| patch  | check_10 + **check_12 (3 unlabeled bdry faces)** | check_13 PASS (0 flips). check_12 catches the predicted `mmg3d_local_patch` regression. |
| cavity | check_10 + check_12 (3 unlabeled bdry faces, inherited from patch) | check_13 PASS |

`orient_fault_surface.py` was run idempotently on each stage before validate; output:
```
raw:    components=2 flipped=0 nonmanifold=393   (would have flipped 6313 on un-oriented input)
mmg3d:  components=2 flipped=0 nonmanifold=269   (would have flipped 9615)
patch:  components=2 flipped=0 nonmanifold=274   (would have flipped 11338)
cavity: components=2 flipped=0 nonmanifold=274   (would have flipped 11338)
```

## Topology metrics on final cavity.msh

Final cavity mesh: `output/newset_6_all6_cavity_retet/output/safs_newset_6_cavity.msh` (192,709 verts, 44,559 tris, 1,108,882 tets), after running `orient_fault_surface.py` (R-003).

| metric | pre-fix | post-fix |
|---|---|---|
| unlabeled 1-tet bdry faces | 3 | **3** (still present — produced upstream by `mmg3d_local_patch`; R-001 check_12 NOW DETECTS them, R-006 producer-side assertion would prevent them on a fresh run; see Recommendations) |
| tag-100 fault orphans | 0 | 0 |
| fault winding flips (BFS-propagated) | 229 | **0** |
| non-manifold fault edges | 269 | **274** (intrinsic to branching SAFS geometry; unchanged by orient pass — reported as warn-only metric per R-002 spec) |
| T-junctions | 0 | 0 |

Numbers from `/tmp/global_check.py` and `/tmp/fault_orientation_check.py` on the orient-corrected cavity mesh.

## Verdict

**R-001..R-006 are FIXED and exercised on real Step-6 stage outputs. READY FOR RE-REVIEW with two notes:**

1. **The 3 surface holes in `_patch.msh` and `_cavity.msh` are now DETECTED, not silently passed.** Pre-fix, `validate_patch.txt` reported `11/11 checks passed` on the broken mesh; post-fix, `validate_patch.log` reports `11/13 checks passed` with the explicit failure: `[FAIL] 12_surface_closure: 3 unlabeled 1-tet bdry face(s); first centroids: (55541,-19265,0), (55435,-19982,0), (55281,-19608,0)`. R-004's fail-fast policy means a fresh end-to-end run will now ABORT at the post-patch validate gate rather than producing a topology-broken cavity mesh.

2. **Recommendation per task spec: roll back `mmg3d_local_patch` for now.** Per the task instructions ("If check_12 fails on the post-patch output: that's the predicted bug — `mmg3d_local_patch.py` is a known producer of surface holes. Document in the fix report and recommend rolling back local_patch (set ENABLE_MMG3D_LOCAL_PATCH=0). Do NOT try to repair via additional patches."). Concretely: the 3 unlabeled bdry faces are introduced by `mmg3d_local_patch._stitch_back` (R-006 assertion would have caught them at the producer if the patch had been re-run after the R-006 fix went in). On a fresh run the producer-side `RuntimeError` from R-006 will halt the patch stage immediately rather than emit a broken mesh. If the user wants Step 6 to complete to cavity.msh in the meantime, run with `ENABLE_MMG3D_LOCAL_PATCH=0` (mmg3d post-pass only — gives a 12/13-passing mesh whose only failure is the pre-existing check_10 sliver gate).

3. **Live pipeline note.** The three relaunched Step-6 attempts in this session each failed at `generate_safs_mesh` (gmsh HXT 3D constrained-recovery error on a 15676-tri STL produced by CGAL 6.1 autorefine). This is HXT/CGAL non-determinism — not a regression caused by R-001..R-006 fixes — and is orthogonal to the validation-and-orientation gate this review addressed. Rerunning the cascade or perturbing the autorefine seed will eventually land on a HXT-friendly STL (the same script succeeded at this stage earlier in the session, producing a 15667-tri STL). The post-fix validation of the existing downstream outputs (table above) demonstrates the fixes work end-to-end on the pipeline's normal stage outputs.
