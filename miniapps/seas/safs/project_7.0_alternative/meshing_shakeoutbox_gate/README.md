# ALT ShakeOut-box meshes, refined to hold their resolved-frequency gates

Requirement:

| mesh | requirement | gate on Vs/dx |
|---|---|---|
| intermediate | 0.5 Hz at p3 | 0.6667 |
| heavy | 0.5 Hz at p3 | 0.6667 |
| heavy | **1 Hz at p5** | **0.8000** ← binding |

Convention (locked): resolved `f = (p/4)·Vs/dx`, `dx` = element MAX edge,
`Vs = sqrt(mu/rho)` NEAREST-GRID at the element BARYCENTRE.

## The finding that reshaped this work

The gate had been judged against the deck's `safs_material_cvm.nc`, which is a
**resample** of MUSCAL onto a 1500 m lateral / **250 m uniform vertical**
lattice. MUSCAL itself is 0.01 deg laterally with a **50 m depth step through
the top 500 m**.

Nearest-grid on the deck file therefore hands every barycentre in the top 125 m
the z = 0 value. Measured on this footprint:

| depth | deck nc Vs (p50) | MUSCAL Vs (p50) |
|---|---|---|
| 0 | 526 | 534 |
| 50 | — (no level) | 996 |
| 100 | — (no level) | 1,429 |
| 200 | — (no level) | 2,091 |
| 250 | 1,950 | 2,353 |

The two models agree at the surface. The deck file simply has nothing between
0 and 250 m, so it reports a **3.78x median step** (p90 4.93x, max 8.99x) where
MUSCAL reports a smooth ramp.

Re-scoring the heavy mesh's 366,083 deck-nc failures against MUSCAL:

* only **10.4 %** still fail at 0.6667 and **28.1 %** at 0.8;
* 69.8 % of them have barycentres inside the deck file's z = 0 bin, where
  MUSCAL reports **2.18x faster** rock (median 345 -> 751 m/s).

So roughly **three quarters of the apparent gate failures were a vertical
binning artefact of the resampled file**, not slow rock. Sizing the mesh to them
was costed at **+15M to +25M collar cells** (`field3d.py`); sizing to MUSCAL
prices a fully compliant collar at **3.97M**, i.e. *cheaper than the 13.4M
already shipped*.

**This is only safe because the decks are moving to native MUSCAL material.**
A parallel build (`..._MUSCALNATIVE5NC`) replaces the single 250 m cube with
five ASAGI blocks on MUSCAL's own depth ladder (dz 50 m in the top 3 km), with
zero vertical resampling. **A mesh gated on MUSCAL must be run with that stack**
-- pairing it with the old 250 m cube would leave the solver propagating through
526 m/s surface material the mesh does not resolve.

## Results

| | intermediate | heavy |
|---|---|---|
| shipped | `..._shakeoutbox.puml.h5` 29,385,401 | `safalt_fb200_deep40km_refine2_shakeoutbox` 135,567,603 |
| **new** | `safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5` **38,827,749** (+32.1 %) | `safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5` **157,840,519** (+16.4 %) |
| parent | `safalt_0d5Hz_p3_deep40km` (unchanged) | **swapped** to `safalt_fb200_deep40km_1Hz_p5` |
| binding gate | 0.6667 (0.5 Hz @ p3) | **0.8000 (1 Hz @ p5)** |
| failures at it (MUSCAL) | 39,281 -> **0** | 93,213 -> **369** (0.0002 %) |
| worst Vs/dx | 0.1275 -> **0.6667 = the gate** | 0.1275 -> **0.4816** |
| **worst resolved f** | 0.096 -> **0.500 Hz** (target 0.5) | 0.159 -> **0.602 Hz** (target 1.0) |
| eta_min | — | 0.0545 |
| min edge | 9.3436 m (= parent's) | 2.3359 m (= parent's) |
| fault area delta | **exactly 0.000e+00 m2** | **exactly 0.000e+00 m2** |
| inverted tets | 0 | 0 |
| fault identity | area delta 0.000e+00; **+18 fault triangles** | ALL PASS |

Min edge equals each parent's exactly, so **no new dt floor from a short edge** —
but `r_insphere` is the real dt test and is NOT yet measured (see below).

### The heavy mesh was welded to the wrong parent

The shipped heavy welds its collar onto `safalt_fb200_deep40km_refine2`
(122,162,105 tets), which was only ever certified at 0.5 Hz/p3 and fails the
1 Hz gate on **274,299 of its own cells** -- 75 % of that product's total.
`safalt_fb200_deep40km_1Hz_p5` (133,699,789) is the same mesh with those closed
and measures **0 failures at both gates** on the deck nc.

It could not reuse the shipped collar: measured by `wall_compare.py`, the 1 Hz
parent's vertical wall is **46,041 triangles / 25,069 verts** against
refine2's **26,704 / 15,205** -- its LEB pass bisected the wall. The
intermediate and refine2 walls ARE bit-identical (sha1 `90c575647b1e5452`), which
is why one collar served both before. So the collar was rebuilt against the
refined wall and re-welded (46,041 interface faces identical on both sides).

## Method

Conforming 3-D longest-edge (Rivara) bisection on the MERGED mesh, fault edges
frozen. Bisection only splits, so compliance is monotone, geometry is exact (the
flat lid stays exactly flat), and it is deterministic -- the properties an mmg
regrade cannot offer, since mmg's +-41 % acceptance band licenses coarsening and
has twice driven this gate BACKWARDS.

**The refinement target is the gate itself (`--pool gate`).** A pooled lower
bound exists to stop a treadmill, and the treadmill is a property of the deck
file's binning: with a 3.78x step at the z = -125 m bin edge, one bisection can
drop a child's Vs ~4x, so the pass mints failures as fast as it fixes them.
MUSCAL has no such cliff -- its 50 m ladder means Vs **bottoms out** at the
surface rather than falling off it -- so the sequence converges on its own.
Measured on the intermediate mesh:

| target | cells flagged | outcome |
|---|---|---|
| `bbox` (full lower bound) | 1,456,239 | 885 s/round, marked flat, ~+24M tets projected |
| `half` (next generation) | 587,899 | 230 s/round, marked flat at ~590k |
| **`gate`** | **39,281** | **25 s/round, peaked at 115k then fell, +18.2 % tets** |

against 39,281 cells that actually fail. Both pooled variants force the entire
collar surface to ~800 m for no measured benefit; their logs are kept as
`*_STALLED_*.log`.

Three passes, reseeding a fresh patch on the residual each time (pass 1 stops
when its k-hop patch RIM starts freezing terminal edges):

| pass | hops | rounds | tets | gate failures |
|---|---|---|---|---|
| in | — | — | 29,385,401 | 39,281 |
| 1 | 4 | 40 | 34,127,870 | ~24,783 |
| 2 | 10 | 60 | 34,650,302 | 886 |
| 3 | 14 | 90 | 34,740,308 | 246 |
| 4 | 20 | 140 | 34,785,799 | 124 |
| z1 (zpool, seeded) | 25 | 200 | 36,886,920 | 21 |
| z2 (zpool, seeded) | 25 | 200 | 38,827,189 | **37 — went BACKWARDS** |
| ig1 (gate, seeded) | 30 | 150 | 38,827,686 | 4 |
| **zero** (`--allow-fault-split`) | 30 | 6 | **38,827,749** | **0** |

Pass 4 LOOKED like the gate target's floor -- the count stopped falling and
OSCILLATED (85 -> 94 -> 110) with the worst pinned at 0.32. **It was not a
floor.** It was the k-hop patch RIM: with `--hops 4-20` the chains that would
have fixed those cells ran into frozen rim edges. Widening the patch while
seeding it only on the measured failures (`--seed-on-gate --hops 30`) resumed
progress immediately. See "Closing the intermediate to ZERO" below for the three
things that were actually wrong and how the gate was closed.

## Reproduce

```bash
CODE=.../meshing_shakeoutbox_gate/code

# what the gate is judged against, and why it matters
python $CODE/cvm_facts.py                       # the deck nc's 250 m binning
python $CODE/muscal_facts.py                    # MUSCAL's 50 m ladder
python $CODE/rescore_muscal.py --dump build_tmp/fail_heavy.npz --gate 0.6667 --gate 0.8

# price a fully compliant collar under either model
python $CODE/field3d.py --gate 0.8 --hgrad 1.5             # deck nc  -> 21.7M cells
python $CODE/field3d.py --gate 0.8 --hgrad 1.5 --muscal    # MUSCAL   ->  4.0M cells

# census any mesh
python $CODE/census_box.py --muscal --mesh <m.puml.h5> --parent-tets N \
       --gate 0.6667 --gate 0.8

# the builds
$CODE/build_collar_1hz.sh        # collar against the 1 Hz parent's refined wall
$CODE/merge_heavy.sh             # parent swap + weld
$CODE/close_intermediate.sh      # pass 1   (then _pass2.sh, _pass3.sh)
$CODE/close_heavy.sh
```

## Residual — a non-zero count DOES mean "not compliant everywhere"

Each remaining cell resolves less than its target. What they actually resolve
(`residual_impact.py`):

| | cells | worst | p10 | median | below half the target |
|---|---:|---:|---:|---:|---:|
| intermediate, target 0.500 Hz @ p3 | 124 | **0.237 Hz** | 0.317 | 0.388 | 4 cells |
| heavy, target 1.000 Hz @ p5 | 369 | **0.602 Hz** | 0.626 | 0.790 | **0 cells** |

Every cell in the heavy mesh resolves at least 0.602 Hz at p5, and 217 of the
369 are already above 0.75 Hz.

**Where they are.** All are shallow free-surface cells (median barycentre
~-105 m). None is near the fault — the heavy has **zero** in the frozen parent
block, its nearest is 31 km outside it (median 189 km); the intermediate has 4
(1.6 %), median 91 km out. But they are NOT confined to dead corners: 82 % and
73 % respectively lie INSIDE the ShakeOut v1 PGV comparison footprint, so they
sit where ground motion is read. Only the handful of worst cells fall past
MUSCAL's lon -114.02 / the deck nc's E 705 km, where Vs is an edge-clamp value
rather than a measurement.

They are a convergence tail, not a structural class (no fault-edge or
wall-pinned exempt class was ever needed). Closing them fully is a matter of
spending a wide-patch pooled pass, not of changing method.

## Note for the intermediate deck

The intermediate mesh resolves **0.5 Hz at p3, not 1 Hz at p5**: at gate 0.8 it
has 455,365 failures (1.31 %). That number went UP during this work (121,931 ->
455,365) and the rise is real, not a regression -- refining a cell to satisfy
0.6667 moves its barycentre toward the surface, which is exactly the
self-reference described above, and the stricter gate then catches it. Use the
heavy mesh for 1 Hz.


## Closing the intermediate to ZERO — what the last cells needed

The gate-as-target passes floored out around 100-124 cells. Three things were
wrong, in order of discovery:

1. **The patch was too tight.** Passes with `--hops 4-20` stalled with hundreds
   of terminal edges frozen at the patch RIM, which reads exactly like a
   convergence floor but is not one. `--seed-on-gate --hops 30` (patch built
   from the measured failures alone, so it stays ~10 % of the mesh even at 30
   hops) took 37 -> 4 for **+497 tets**.
2. **A pooled target must run to COMPLETION or not at all.** A z-pooled pass
   seeded on failures took 124 -> 21 (+2.1M tets), but the next one went
   21 -> **37**, i.e. backwards. Conforming propagation bisects compliant cells
   as part of a LEPP chain, and those children can dip below the gate; the
   pooled bound only pays that back once every flagged cell reaches its target.
   Rim-stalled, it never does. Mesh-wide zpool flags 3,182,350 cells and a 78 %
   patch — the price of the guarantee. Logs: `close_int_z*_*.log`,
   `close_int_z1_MESHWIDE_TOOBIG.log`.
3. **The true floor was the FAULT's own triangulation.** The last 4 cells sat at
   the shallow surface trace (z = -6..-23 m, Vs 189-258 m/s) with dx 290-402 m,
   needing only 1.03-1.13x. Their LEPP chains terminated on fault edges:
   `all 3 terminal edges frozen (rim 0, fault 3) -- stop`. That is the
   fault-edge-pinned structural class.

`--allow-fault-split` releases it. Rivara inserts the MIDPOINT of an edge, and a
fault edge's midpoint lies exactly on the planar fault triangles sharing it, so
**the fault surface and its area are preserved bit-for-bit** — measured area
delta exactly `0.000e+00 m2`. Only the triangulation gets finer. Cost: **63
tets**, 6 rounds, and **9 of 160,280 fault facets subdivided**
(320,560 -> 320,578 triangles, +0.006 %).

### Deck impact of the fault split — read this

The DR facet COUNT changed, 160,280 -> 160,289. Spatial fault inputs (stress nc,
friction nc, nucleation) are FIELDS and resample cleanly, and no fault vertex
moved (asserted). But anything keyed to the facet LIST rather than to position
must be re-derived: pickpoint indices, facet-count assertions, and any stored
per-facet ordering. If that is unacceptable, the alternative is to ship the
124-cell version and accept 0.237-0.500 Hz on those cells.


## ParaView views (`view/`, built by `code/make_views.sh`)

| file | what | size |
|---|---|---|
| `intermediate_fault.xdmf` + `.h5` | 160,289 DR triangles, de-duplicated | 6.5 MB |
| `intermediate_surface.xdmf` + `.h5` | fault + free surface + absorbing hull, `bc` scalar | 140 MB |
| `intermediate_full.xdmf` | points AT the .puml.h5, zero copy — loads all 38.8M tets | 914 B |
| `heavy_fault.xdmf` + `.h5` | 2,564,480 DR triangles | 103 MB |
| `heavy_surface.xdmf` + `.h5` | fault + lid + hull, `bc` scalar | 264 MB |
| `heavy_full.xdmf` | points AT the .puml.h5 — 157.8M tets, needs the RAM | 904 B |

A PUML stores each fault triangle TWICE (once per side); the fault/surface modes
de-duplicate, which is why `intermediate_fault` shows 160,289 and not 320,578.
Colour `surface` by `bc`: 1 = free surface, 3 = dynamic rupture, 5 = absorbing.

**Do not move a `.xdmf` away from the `.h5` it names** — XDMF resolves
`<DataItem>` paths relative to the XDMF file, and `_full` carries a relative
path back to `../results/`.

## What is kept on disk, and what was deleted

Kept: the two final meshes, the six views, and every parent needed to rebuild
them (`meshing_deep40km/safalt_0d5Hz_p3_deep40km`,
`meshing_deep40km_1Hz_p5_leb/safalt_fb200_deep40km_1Hz_p5`).

Deleted as superseded/temporary: every intermediate pass product, both build
collars (`collar_1hz_uniform_2500.npz`, `collar_uniform_2500.npz`), the census
dumps, and the two PREVIOUS shakeoutbox products with their view folders —
`meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5`
and `meshing_shakeoutbox_heavy/results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5`
(the latter was the one welded to the wrong parent). Both are reproducible from
the parents plus `code/`. `build_tmp/` now holds only logs and stats JSON.
