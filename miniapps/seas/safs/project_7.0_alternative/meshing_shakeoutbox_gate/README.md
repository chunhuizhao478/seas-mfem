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
| **new** | `safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5` **34,785,799** (+18.4 %) | `safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5` **157,840,519** (+16.4 %) |
| parent | `safalt_0d5Hz_p3_deep40km` (unchanged) | **swapped** to `safalt_fb200_deep40km_1Hz_p5` |
| binding gate | 0.6667 (0.5 Hz @ p3) | **0.8000 (1 Hz @ p5)** |
| failures at it (MUSCAL) | 39,281 -> **124** (0.0004 %) | 93,213 -> **369** (0.0002 %) |
| worst Vs/dx | 0.1275 -> **0.3158** | 0.1275 -> **0.4816** |
| **worst resolved f** | 0.096 -> **0.237 Hz** (target 0.5) | 0.159 -> **0.602 Hz** (target 1.0) |
| eta_min | — | 0.0545 |
| min edge | 9.3436 m (= parent's) | 2.3359 m (= parent's) |
| fault area delta | **exactly 0.000e+00 m2** | **exactly 0.000e+00 m2** |
| inverted tets | 0 | 0 |
| fault identity | ALL PASS | ALL PASS |

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
| 4 | 20 | 140 | **34,785,799** | **124** |

Pass 4 is where the gate target reaches its floor: the count stopped falling and
began to OSCILLATE (85 -> 94 -> 110) with the worst pinned at 0.32. That is the
treadmill appearing at the extreme tail, in the few columns slow enough that a
bisected child lands in slower material than its parent.

A 5th pass was attempted with a HYBRID -- seed the patch on the ~100 measured
failures (so it is 0.30 % of the mesh instead of 55-66 %) but drive it with the
pooled bound, which is the treadmill's cure (`--seed-on-gate --pool half`). It
did not help: the tight patch froze ~375 terminal edges at its own rim and the
pooled count went sideways. Log kept as `close_tail_int_RIMSTALLED.log`. Closing
the last ~100 cells needs a wide patch AND the pooled target together, which is
the expensive combination this campaign exists to avoid.

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
