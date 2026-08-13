# ALT ShakeOut-D heavy mesh — 1 Hz at p5

**Product:** `results/safalt_shakeoutD_heavy_1Hz_p5.puml.h5` (3.88 GB)

```
87,240,911 tets   16,149,797 verts
volume 12,480,379 km3  -- EXACT match to the ShakeOut-D domain
inverted / zero-volume tets: 0
```

## Acceptance

| gate | result |
|---|---|
| 1 Hz at p5 (`Vs/dx >= 0.8`, native MUSCAL, dx = MAX edge) | **35,353 of 87,240,911 below gate = 0.041 %** |
| worst cell | `Vs/dx` 0.2990 -> resolved **0.3737 Hz** |
| failing cells | 1st pct 0.604 Hz, median 0.908, 99th 0.998 — all marginal |
| median resolved | 13.4 Hz |

PUML structure — every check PASS:

```
P1 BC3 fault        2,564,469 distinct, ALL interior (x2), 0 exterior
   both incident tets tagged: 2,564,469 / 2,564,469
P2 BC1 free-surface 3,511,801 all exterior
P2 BC5 absorbing       78,177 all exterior
P3 untagged exposed faces (holes): 0
P4 inverted/zero-volume tets: 0
```

## Fault preserved verbatim

| | value |
|---|---|
| facets | 2,564,469 (parent surface 2,564,480; 11 lost = pinholes) |
| area | **6,537.69 km2** vs 6,537.70 extracted — 1.5e-6 relative |
| max vertex displacement through mmg | 2.046e-11 m (tol 1e-3) |
| centroid match after mmg | 9.601e-10 m |
| lost area | 5,307 m2 = **8.1e-05 %** |

## Build path

| stage | result |
|---|---|
| PLC | 4,603,808 facets, Euler V−E+F = 2, all 8,137 trace edges embedded |
| base fill (`tetgen -pY`) | 9,613,853 tets, volume exact, 6 pinholes |
| mmg (`-hgrad 1.3 -hausd 30 -hmin 115 -hmax 5000`) | 63,044,732 tets, 8h47m |
| rv repair (20 iters) | rv>2 559,382 -> 388,640; rv>3 261,973 -> 219,733; med 1.4209 -> 1.3963 |
| gate close (LEB, 40 rounds) | 414,849 -> 35,353 below gate (**-91.5 %**), +24.2 M tets |

## Known deviations — read before using

1. **Gradation is g ≈ 0.30, not the specified 0.15.** mmg ran at `-hgrad 1.3`,
   which is the ratio `1+g`. Reaching 0.15 needs a second mmg pass against
   `build_tmp/target_g015.npz` (built, pair-tested, delivered gradation
   <= 0.1499 against a 0.15 spec) at ~9-12 h and ~27 GB. **Deliberate user
   decision** after the LEB route to g=0.15 was shown to diverge.
2. **Gate is 99.96 % closed, not 100 %.** 35,353 cells remain, all marginal
   (median 0.908 Hz). The LEB pass hit its 40-round cap while still improving
   ~200 cells/round; growth had nearly stopped (+69 k tets over the last 2
   rounds), so continuing is cheap but slow.
3. **Slivers: 1.2 M below eta 0.1, 100 % fault-adjacent** (97.4 % within 250 m,
   zero beyond 2 km). rv repair reached what it could — 18 of 20 iterations took
   no step, i.e. it had converged. Structural: tets with a face on the frozen
   fault have 3 of 4 vertices immovable. ALT's fault carries 352 needle triangles
   with a 2.34 m minimum edge (a sibling mesh's is 21.14 m), and no well-shaped
   tet exists on a 2.34 x 98 x 99 m triangle.
   **Context: the deployed small-domain mesh carrying this same fault reached
   eta_min 0.0307 with 8 cells below 0.05, via gmsh+LEB rather than tetgen CDT +
   mmg.** The sliver population is a property of the base-fill route, not the
   geometry. Worth revisiting only if dt proves unacceptable.

## Two measurement traps hit during this build

* **`make_view_xdmf.py` without `--muscal` silently scores on the DECK cube.**
  `Material` returns early when `muscal_nc is None`, leaving `self.M` unset, and
  `at()` falls back to the deck despite `source="muscal"`. The deck's 250 m
  z-binning manufactures shallow failures: it reported **942,375** cells below
  gate where the true MUSCAL count is **35,353** — a 27x overstatement that
  briefly looked like the gate pass had failed. Now raises instead.
* **`leb_gate_close.py`'s per-round `gate-fail` counts only the ACTIVE set.**
  `fail` is zero-initialised and filled only for `idxs`; after round 0
  `active = fail.copy()`. Round 0 is a whole-mesh count (414,849, independently
  reproduced); later rounds are not. Verify the product independently — here the
  final whole-mesh count (35,353) happened to agree with the last active-set
  count (35,501), but that is not guaranteed.

## Reproduce

```
code/run_mmg_s1.sh                    # base -> 63.0 M   (--keep-medit MANDATORY)
code/medit_to_puml.py build_tmp/s1.mmg_out.mesh build_tmp/..._s1.puml.h5
code/run_finish.sh 20                 # rv repair + gate close, chained
code/check_puml_faces.py results/safalt_shakeoutD_heavy_1Hz_p5.puml.h5
code/make_view_xdmf.py --mesh ... --order 6 --cvm ... --muscal ...   # --muscal REQUIRED
code/quicklook.py --view view --tag alt_final
```

Views: `view/alt_final_{full,surface,fault}.xdmf` + `alt_final_quicklook.png`.
