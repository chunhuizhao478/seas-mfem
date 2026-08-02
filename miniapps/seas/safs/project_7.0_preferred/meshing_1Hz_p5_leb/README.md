# meshing_1Hz_p5_leb (PREFERRED)

Refines `safpref_fb200_refine2_gatefix.puml.h5` so it resolves **1 Hz at p5
everywhere**, by locating the failing cells and applying **local conforming Rivara
longest-edge bisection** — skill architecture 4, the same recipe as the ALT build
in `project_7.0_alternative/meshing_deep40km_1Hz_p5_leb/`.

## Product

    results/safpref_fb200_refine2_gatefix_1Hz_p5.puml.h5
      124,232,100 tets / 21,504,494 verts
      5,485,395,952 bytes
      md5 963e0a2e2669543f5b94b46154383a10

    resolved everywhere:  1.0000 Hz @ p5   |   0.6000 Hz @ p3
    gate failures:        0  (parent: 428,060)
    run cost:             +2.7 % over the parent (dt_min unchanged)

Full numbers and the deck verification: `QUALITY_REPORT_1Hz_p5.md`.

## The one thing that differs from the ALT build

**The PREFERRED lid carries OFF-DATUM vertices** where the ALT lid is exactly flat.
It is *not* topographic — measured 2026-08-01, only 2,961 of 975,156 top vertices have
|z| > 1e-6 and 2,417 of those are sub-0.1 m noise; the large ones are the **fault trace**
where the fault daylights (73 are fault vertices, same 10 coordinates on the small and
intermediate meshes). Flattening them would move the fault, so it was declined.
Two checks still had to change, and both are now in this folder's `code/`:

* `check_fault_identity.py` — F4 can no longer assert "exactly flat". It now checks
  the free-surface **total area** and **z-range** are unchanged. LEB inserts edge
  midpoints, and the midpoint of an edge of a planar triangle lies ON that triangle,
  so a correct refinement leaves the lid area bit-unchanged. Measured: agrees to
  **1.5e-16 relative**. (On a flat lid this reduces to the same statement, so the
  check now serves both lineages — prefer this version.)
* `check_receivers_local_top.py` — **new**. Receiver containment must be barycentric
  under the **local** top triangle, not a 2-D footprint test, because SeisSol v1.1.3
  silently DROPS receivers above their local free surface. Vertex proximity is not
  enough either: the far-field lid is kilometres coarse, so the nearest top vertex
  can belong to a triangle that does not contain the point. Result: 93,700/93,700
  below the local lid, min clearance 0.290 m — which **reproduces the parent's
  documented 0.290 m exactly**, independently confirming the lid did not move.

## Also worth knowing

* The census found **9 fault-vertex-pinned** cells (ALT had 0). Under this method
  that is **not** a structural floor: the refiner freezes fault EDGES, not fault
  VERTICES, and bisecting a non-fault edge that merely touches a fault vertex does
  not move the fault. All 9 were refined; the gate closed to zero.
* ⚠ **Node sizing at p5 needs resolving before any run.** A circulating estimate
  uses "ORDER 5 = 35 basis functions = 1.75×" as the p5 factor, but 35 basis
  functions is degree **4**. A true p5 run is ORDER 6 = 56 basis functions =
  **2.80×** ORDER 4. See §6 of the quality report.

## Layout

    code/census_1Hz_p5.py             locate + classify gate failures (barycenter rule)
    code/leb_refine_1Hz.py            the refiner (chunked, vectorized; fault edges frozen)
    code/check_fault_identity.py      F1-F5, lid area+range check (any lid)
    code/check_receivers_local_top.py Stage F receivers under the LOCAL lid
    code/puml_io.py                   PUML read/write + packed BC word
    results/                          product + census/stats json + mesh.md5
    build_tmp/                        run logs

## Reproduce / knobs

See §8 of `QUALITY_REPORT_1Hz_p5.md`. Two passes are needed: pass 1 (`--hops 3`)
stops when its remaining terminal edges hit the patch rim, leaving 9,548 gate
failures; pass 2 (`--hops 8`) reseeds on just those and closes to zero.
`--gate` 0.8 = 1 Hz @ p5, 0.6667 = 0.5 Hz @ p3.

## Memory

Peak RSS 16.6 GB (pass 1, 96 min), 12.1 GB (pass 2, 14 min). Do **not** run a second
mesh-refine job of this size concurrently — census stages (~4 GB, 25 s) are safe to
overlap, refine stages are not.

## Consuming deck

`~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_PREFERRED_THERMAL_CASE2_fb200ref2gatefix_plast_phi30_40_fw15pregate_nw08_k1p75_mw783_attenuation`
currently points at the parent (`mesh_preferred_fb200_ref2_gatefix.puml.h5`). All
fault-referenced inputs transfer unchanged (fault area delta exactly 0, DR facets
identical, 93,700/93,700 receivers contained, extent identical). To switch it, update
`MeshFile` plus `MESH`/`MESH_BYTES`/`MESH_MD5` in the sbatch — but resolve the p5
node-sizing question first.
