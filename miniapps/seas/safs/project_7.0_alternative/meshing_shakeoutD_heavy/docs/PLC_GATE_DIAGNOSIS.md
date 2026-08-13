# ALT ShakeOut-D PLC gate: `recoversubfaces` internal error is a `-d` artifact

`isolate_selfint.py` on the assembled ALT PLC returned

```
The input surface mesh is correct.
  fault only              2,564,480 facets   (Warning: 2 segments are not recovered)
The input surface mesh is correct.
  hull only               2,039,328 facets
  top only                1,830,914 facets  ->  All vertices are coplanar
  hull + fault            4,603,808 facets  ->  Internal TetGen error within `recoversubfaces`
```

**This does not block the build.** Two things had to be separated first.

## 1. "Failed to tetrahedralize" is the PASS signature, not a failure

`isolate_selfint.py` runs bare `switches="d"`. `-d` only *detects*; it never emits a
tetrahedralization, so the wrapper always raises. PREFERRED's **final, known-good**
PLC -- the one that went on to fill -- logged exactly the same pair of lines
(`project_7.0_preferred/.../build_tmp/selfint7.log`):

```
The input surface mesh is correct.
  hull + fault            4,799,610 facets  ->  Failed to tetrahedralize...
```

So ALT's `fault only` and `hull only` results are **clean**. Only `hull + fault`
differs, and it differs by crashing *before* reaching a verdict.

## 2. Every structural hypothesis was tested and refuted

| hypothesis | test | result |
|---|---|---|
| bad fault/hull weld, near-duplicate vertices | KD weld audit | 8,137 welded at **0.000e+00 m**, all distinct; **0** vertex pairs < 0.5 m in the whole PLC |
| non-manifold (branch) edges break recovery | edge valence census | ALT **840**; PREFERRED **1,048** and it filled -> refuted |
| T-junction (vertex on a non-incident edge) | exhaustive local vertex-vs-edge scan around all needles + junctions | **0** within 1 m |
| needle/sliver triangles (min edge 2.34 m vs PREFERRED's 21.1 m) | see section 3 | recovered fine -> refuted |
| shallow "flat wedge" where the fault leaves z=0 | count tris with >=2 verts at z=0 and abs(nz)>0.5 | ALT **232** (max 0.72); PREFERRED **2,340** (max 0.98) and it filled -> refuted |
| fault poking above the free surface | z audit of unwelded fault verts | max z = **-5.076 m**, 0 above -1e-9 |
| trace edges in the z=0 plane never embedded (**PREFERRED's actual fatal defect**) | every fault edge with both ends welded vs top-surface edge set | **0 of 8,137 missing** |
| fault touching a wall/bottom | rotated-frame containment | clearances 59.8-165.9 km lateral, 60.9 km bottom |

## 3. The decisive test: reduced `-pY` on the worst geometry

The full fill is the one stage that must not overlap another heavy job, so the real
switches were run on small closed boxes around the worst geometry instead --
identical code path, ~1 % of the memory (`tmp/localbox.py`).

| box | fault tris | min edge | min q | needles | non-manifold edges | `-pY` result |
|---|---|---|---|---|---|---|
| sliver cluster, deep | 17,772 | 2.590 m | 0.048 | 36 | 0 | **OK** 40,507 tets, 10 Steiner |
| sliver cluster, to surface | 9,553 | 2.590 m | 0.048 | 36 | 0 | **OK** 23,442 tets, 9 Steiner |
| junction comp 0 | 34,974 | 5.011 m | 0.120 | 0 | 130 | **OK** 92,616 tets, 10 Steiner |
| junction comp 2 | 195,885 | 5.216 m | 0.078 | 36 | **403** | **OK** 503,065 tets, 52 Steiner |

533 of the 840 non-manifold edges (63 %) and the worst slivers all recover under
`-pY`. The Steiner points are **interior**; `-Y` forbids them only on input facets,
so the frozen fault triangulation is preserved in every case.

**Conclusion:** the `recoversubfaces` crash is a robustness failure of tetgen's `-d`
diagnostic path on a 4.6 M-facet PLC, not a defect in the geometry. Proceed to
Stage 2 (`fill_domain.py --minratio 0 --hmax 0`), which is the real gate.

## Note on ZTOP = -60

Worth recording because it looked like the culprit and is not. The deployed ALT mesh
has 7,785 fault-boundary vertices that already daylight plus a 390-vertex buried
shallow tip (-5.03 .. -142.46 m); exactly **352** of those sit above -60, which is
precisely what `extract_fault_surface.py` snapped (8,137 - 7,785 = 352). So ZTOP=-60
does daylight part of the tip the original build left blind. But the resulting flat
wedges are 10x *fewer* and less extreme than PREFERRED's, which filled without
trouble, so this is not a build blocker and ZTOP=-60 stands as chosen.
