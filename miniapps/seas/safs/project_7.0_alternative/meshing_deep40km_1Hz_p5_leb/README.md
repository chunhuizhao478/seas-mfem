# meshing_deep40km_1Hz_p5_leb

Refines `safalt_fb200_deep40km_refine2.puml.h5` (the fault-band L2 mesh) so it
resolves **1 Hz at p5 everywhere**, by locating the cells that fail the gate and
applying **local conforming Rivara longest-edge bisection** — not a global remesh.

## Product

    results/safalt_fb200_deep40km_1Hz_p5.puml.h5
      133,699,789 tets / 23,123,324 verts
      5,902,955,432 bytes
      md5 554ab7448c6a2d02388f01686e333279

    resolved everywhere:  1.0000 Hz @ p5   |   0.6000 Hz @ p3
    gate failures:        0  (parent: 274,299)
    run cost:             +1.7% over the parent (dt_min unchanged)

Full numbers, regressions and the failed approaches: `QUALITY_REPORT_1Hz_p5.md`.

## Why bisection

Two earlier mmg size-map regrades of this mesh **regressed** the gate
(0.6667 → 0.4104 / 0.4137): mmg's ±41% acceptance band licenses coarsening, so
compliant cells drift over tolerance and the pass mints failures as fast as it
fixes them. Bisection only ever splits, inserts midpoints only (so geometry is
exact and the flat top stays exactly flat), is conforming by construction, and is
deterministic. It also left every worst-case shape metric bit-identical to the
parent.

## The one non-obvious design decision

The refinement **target** is not the gate. The gate is `Vs` nearest-grid at the
barycenter (locked convention). Targeting that directly *treadmills*: the CVM is
nearest-grid on a 250 m vertical lattice with a ~4× Vs step at z = −125 m, so
bisecting a straddling cell throws one child into the slow bin needing 4× smaller
dx. Measured: stalled at ~65,000 failures, −1.8%/round while +157,000 tets/round.

The shipped target is `min Vs over the cell's own VERTICAL extent`, sampled at the
barycenter's (x,y). It is a lower bound on any child's measured Vs, hence
refinement-stable. Pooling in **z only** matters: a full 3-D pool cost +55.8% tets
on the small ALT mesh (it takes the slowest Vs anywhere in a multi-km cell), while
the CVM is 1500 m laterally vs 250 m vertically, so bin changes are essentially
always vertical.

## Layout

    code/census_1Hz_p5.py        locate + classify gate failures (the GATE, barycenter rule)
    code/leb_refine_1Hz.py       the refiner (chunked, vectorized LEB; fault edges frozen)
    code/check_fault_identity.py F1-F5: fault multiset, area, BC round-trip, flat top, inverted
    code/puml_io.py              PUML read/write + packed BC word (copied from meshing_deep40km)
    results/                     product + census/stats json + quality_compare.json
    build_tmp/                   run logs, including the two STALLED runs kept as evidence

## Reproduce

```bash
CVM=<deck>/safs_material_cvm.nc
P=<...>/meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5

python code/census_1Hz_p5.py  --mesh $P --cvm $CVM --gate 0.8 --out results/census_round0
python code/leb_refine_1Hz.py --mesh $P --cvm $CVM --out pass1.puml.h5 \
       --gate 0.8 --hops 3 --max-rounds 80  --stats results/leb_stats.json
python code/leb_refine_1Hz.py --mesh pass1.puml.h5 --cvm $CVM \
       --out results/safalt_fb200_deep40km_1Hz_p5.puml.h5 \
       --gate 0.8 --hops 8 --max-rounds 120 --stats results/leb_stats_p2.json
python code/census_1Hz_p5.py        --mesh results/...1Hz_p5.puml.h5 --cvm $CVM --gate 0.8 --out results/census_final
python code/check_fault_identity.py --parent $P --new results/...1Hz_p5.puml.h5
```

Two passes are needed: pass 1 terminates when its remaining terminal edges hit the
patch rim (`hops=3`), leaving 2,253 gate failures; pass 2 reseeds a fresh patch on
just those and runs with `froze rim 0` throughout.

**Knobs.** `--gate` 0.8 = 1 Hz @ p5, 0.6667 = 0.5 Hz @ p3. `--hops` = patch radius;
raise it if the log reports terminal edges frozen at the rim. `--max-rounds`.

## Memory

Peak RSS 20.0 GB (pass 1), 11.8 GB (pass 2) on a 38.7 GB machine. Every geometry
pass is chunked — an unchunked `P[T[:, EI]]` on 122M tets is a 17.6 GB temporary and
will swap-thrash the machine (it did, once). Do **not** run a second mesh job of
this size concurrently; the census stages (~4 GB, 30 s) are safe to overlap, the
refine stages are not.

## Consuming deck

`~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_fbrefine2_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km`
currently points at the **parent** refine2 mesh. Every fault-referenced input
transfers unchanged (fault area delta exactly 0, DR facets identical, 97,711/97,711
receivers contained, extent identical). To switch it, update `MeshFile`,
`MESH`/`MESH_BYTES`/`MESH_MD5` in the sbatch, and re-derive nothing else — but note
the on-fault output and node sizing are unchanged because DR facets and dt_min are
identical.
