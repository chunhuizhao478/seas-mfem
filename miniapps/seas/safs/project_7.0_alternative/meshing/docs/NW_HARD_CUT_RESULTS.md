# NW Hard-Cut Results — ALT6 Long-Strip Fault

**Date run:** 2026-05-08
**Driver:** `code_preprocess/nw_cut_strip.py --batch --out-dir data_cutnwfault -v`
**Reference (defines NW cut anchor):**
`project_7.0_preferred/data_cleanfreesurf/SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl`
**Inputs (read-only, MD5-verified unchanged after the run):**
`project_7.0_alternative/data_cleanfreesurf/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_{500,1000,2000}m_clean_clip.stl`
**Outputs:**
`project_7.0_alternative/data_cutnwfault/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_{500,1000,2000}m_clean_clip_nwcut.stl`

## Cutting plane (constant across all three resolutions)

| Quantity | Value |
| --- | --- |
| NW direction (horizontal unit normal) `n_h` | `(-1, +1, 0) / √2 = (-0.70710678, 0.70710678, 0)` |
| Anchor (preferred mesh's NW-most vertex) | `x = 365072.300 m, y = 3838982.000 m, z = 0.000 m` (UTM 11 N) |
| Anchor NW projection `s_anchor = c` | `2 456 425.1061 m` |
| Half-space kept | `n_h . p ≤ c`, i.e. `(-x + y) ≤ -x_anchor + y_anchor` |
| Pre-cut alt strip max NW projection | `2 692 628.682 m` (≈ 236 km of NW overshoot removed) |

The cutting plane is purely vertical (parallel to the Z axis); only the
horizontal `(x, y)` projection is used for classification.

## Per-resolution end points and counts

For each output mesh:
- **NW endpoint** = the cleaned mesh's vertex with the largest NW
  projection `s = (-x + y) / √2`. By construction this lies on (or
  immediately SE of, within ~1 nm) the cutting plane.
- **SE endpoint** = the cleaned mesh's vertex with the smallest NW
  projection. This is the far south-east corner of the surviving strip.
- **NW–SE strip length** = `s_NW − s_SE` (a 1-D length along the NW
  direction in the horizontal plane).

### 500 m resolution

| | x [m UTM 11 N] | y [m UTM 11 N] | z [m, depth] | s = (-x + y)/√2 [m] |
|---|---:|---:|---:|---:|
| **NW endpoint** | `364978.862` | `3838888.562` | `-8074.304` | `2 456 425.106` |
| **SE endpoint** | `618202.100` | `3692278.000` | `-54.000`   | `2 173 699.915` |
| **NW–SE strip length** | | | | **`282 725.2 m`** (≈ 283 km) |

- Triangles in / out: `71 646 → 42 389` (`29 257` removed, ~40.8 %)
- Cleaned vertex count: see file (output written with `%.6f` precision)
- Drift snap: 63 vertices within 1 m of the plane projected onto it
- Post-cut bbox: `x ∈ [364865.6, 622329.5], y ∈ [3692278.0, 3838985.0], z ∈ [-16607.6, 0.0]`
- File: `data_cutnwfault/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m_clean_clip_nwcut.stl` (11.03 MB)

### 1000 m resolution

| | x [m UTM 11 N] | y [m UTM 11 N] | z [m, depth] | s = (-x + y)/√2 [m] |
|---|---:|---:|---:|---:|
| **NW endpoint** | `364875.015` | `3838784.715` | `-11429.108` | `2 456 425.106` |
| **SE endpoint** | `618202.100` | `3692278.000` | `0.000`       | `2 173 699.915` |
| **NW–SE strip length** | | | | **`282 725.2 m`** (≈ 283 km) |

- Triangles in / out: `18 056 → 10 620` (`7 436` removed, ~41.2 %)
- Drift snap: 32 vertices within 1 m of the plane projected onto it
- Post-cut bbox: `x ∈ [364830.8, 622329.5], y ∈ [3692278.0, 3838944.9], z ∈ [-16607.4, 0.0]`
- File: `data_cutnwfault/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_1000m_clean_clip_nwcut.stl` (2.76 MB)

### 2000 m resolution

| | x [m UTM 11 N] | y [m UTM 11 N] | z [m, depth] | s = (-x + y)/√2 [m] |
|---|---:|---:|---:|---:|
| **NW endpoint** | `364915.675` | `3838825.375` | `-2480.578` | `2 456 425.106` |
| **SE endpoint** | `618202.100` | `3692278.000` | `0.000`     | `2 173 699.915` |
| **NW–SE strip length** | | | | **`282 725.2 m`** (≈ 283 km) |

- Triangles in / out: `4 574 → 2 697` (`1 877` removed, ~41.0 %)
- Drift snap: 17 vertices within 1 m of the plane projected onto it
- Post-cut bbox: `x ∈ [364748.6, 622329.5], y ∈ [3692278.0, 3838955.1], z ∈ [-16607.5, 0.0]`
- File: `data_cutnwfault/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m_clean_clip_nwcut.stl` (0.70 MB)

## Cross-resolution observations

- All three resolutions report the **same SE endpoint** (the original SE
  corner of the alt strip is unaffected by the NW cut). The horizontal
  coordinates `(618202.100, 3692278.000)` are common across all three.
- All three NW endpoints have the same NW projection `s = 2 456 425.106 m`
  (= the plane constant `c` to within `5e-10 m`), but their `(x, y, z)`
  differ between resolutions because each mesh interpolates a slightly
  different point onto the plane (the cut creates new vertices at the
  intersection of the plane with whichever triangle straddles it; the
  three meshes have different straddling triangles at this resolution).
- The NW endpoint's `z` is **not** at the surface: the NW-most vertex of
  each cut mesh sits at `z ≈ -2.5 km` (2000 m), `-8.1 km` (500 m), or
  `-11.4 km` (1000 m) — i.e. on the down-dip edge of the cut face. This
  is geometrically correct: the cutting plane is vertical, so its
  intersection with a dipping fault is a vertical line, and the
  surviving NW edge of the fault is a 1-D segment from `z ≈ 0` down to
  the deepest point of the fault that satisfies `s = c`.
- The strip length along the NW–SE axis is invariant across
  resolutions: `~282 725.2 m` (about 283 km).
- ~41 % of input triangles were removed in every resolution
  (consistent with each resolution's mesh covering the same physical
  surface).

## Validation

- All three input STLs have unchanged MD5 sums after the run:
  ```
  4846dd687e24062b401f37bafad27bbd  …ALT6_1000m_clean_clip.stl
  63995bd34a5e781749e930f7dede382d  …ALT6_2000m_clean_clip.stl
  198f4e2ff0003df05802db386e0a10bb  …ALT6_500m_clean_clip.stl
  ```
- Post-cut `s_max ≤ 2 456 425.1061 + 5e-10 m` for every output ⇒ no
  surviving vertex sits NW of the preferred mesh's NW corner, to within
  float64 cancellation noise.
- All three outputs are ASCII STL (first line `solid nwcut`).
- 28 / 28 unit tests pass (`pytest -q` in
  `project_7.0_alternative/code_preprocess/`).

## Reproducing

```bash
conda activate pythonenv
cd /Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/code_preprocess
python nw_cut_strip.py --batch \
    --out-dir /Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_alternative/data_cutnwfault \
    -v
```

Inspect the anchor without cutting:

```bash
python nw_cut_strip.py --print-anchor /Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/project_7.0_preferred/data_cleanfreesurf/SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl
# → anchor: x=3.650723e+05  y=3.838982e+06  z=0.000000e+00  s=2.456425e+06
```
