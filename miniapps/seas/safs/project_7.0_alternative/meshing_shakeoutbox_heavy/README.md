# ALT heavy mesh, footprint extended to the ShakeOut 2008 box

`results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5`
— 148,966,177 tets / 26,118,137 verts, 6,585,486,464 B.

Parent: `meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5`
(122,162,105 tets), the mesh behind SeisSol run `7885821`.

| | parent | new |
|---|---|---|
| footprint | 253.2 x 420.4 km rotated 30 deg, **106,476 km2** | axis-aligned **706.619 x 464.225 km = 328,030 km2** |
| depth | -40,000 m | unchanged |
| tets | 122,162,105 | 148,966,177 (+26,804,072, **+21.9 %**) |
| fault | 2,564,480 triangles | **identical** |
| free surface | still exactly flat at z = 0 | 5,220,982 facets |

## The collar is SHARED with the intermediate mesh

`meshing_shakeoutbox_intermediate/build_tmp/collar_deep_v6.npz` is welded onto
BOTH deep parents. That is sound because they carry a **bit-identical side
wall**: 26,704 vertical absorbing triangles / 15,205 vertices, verified as an
exact coordinate multiset. `refine2` only ever touched the fault band, so the
far field it inherited from `safalt_0d5Hz_p3_deep40km` is untouched — measured
independently as 10,690,093 vs 10,695,022 cells in the >40 km-from-fault
columns, with an identical per-Vs-bin table.

`merge_collar.py` relocates the wall on whichever parent it is handed, by exact
coordinate key, never by the vertex or tet indices stored in the collar npz —
the two parents number their vertices completely differently (20.9M vs 3.2M).

Build the collar once (see `../meshing_shakeoutbox_intermediate/README.md` for
the sizing sweep and why v6 was chosen), then:

```bash
BASE=$PWD/../..
python code/merge_collar.py \
  --parent $BASE/meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5 \
  --collar $BASE/meshing_shakeoutbox_intermediate/build_tmp/collar_deep_v6.npz \
  --out results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5
```

## Weld and BC census (from the merge)

```
WELD OK: 26,704 interface faces identical on both sides
BC census: fault 5,128,960 (=2x2,564,480), free surface 5,220,982, absorbing 74,818
hull identity: tagged 10,424,760 == hull 5,295,800 + fault 5,128,960 -> OK
collar signed volume: 0 negative
```

`verify_extension.py` is run with `--skip-hull` on this tier: the independent
hull census sorts 4x nt faces, which is ~12 GB *on top of* both meshes at 149M
tets. The O(1) identity above is the same statement and was cross-validated
against the full census on the small and intermediate meshes (177,106 and
5,191,996 hull faces, both exact).

## Known limitation carried over from the collar

The collar holds the 0.5 Hz p3 gate everywhere except a thin near-surface
layer: **363,608 of its 26,804,072 cells (1.357 %) are below Vs/dx = 0.6667**,
89.7 % of them in the top 500 m, worst 0.3163 (0.237 Hz at p3). Below 500 m the
collar is essentially clean. The parent's own 122M cells are bit-identical and
keep their compliance. Collar min edge is 16.08 m against the parent's 9.34 m,
so **the timestep floor is unchanged**.

See `../meshing_shakeoutbox_intermediate/README.md` for the full sizing sweep,
why longest-edge bisection could not close this, and the mmg route that was not
attempted.

## Cost

+21.9 % tets. Quote the LTS cost, not the tet count: the added cells are
shallow far-field and join COARSE LTS clusters, and `dt_min` is set by the
fault band, which is untouched. The parent's DR facet count and its minimum
insphere are unchanged, so the clustered rate-2 LTS cost should rise by far
less than 21.9 % — **measure it before quoting a node-hour figure.**
