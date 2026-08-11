# QUALITY REPORT — ALT ShakeOut-box meshes refined to their frequency gates

Scope: ALT **intermediate** (0.5 Hz @ p3) and ALT **heavy** (0.5 Hz @ p3 AND
1 Hz @ p5). Convention locked: resolved `f = (p/4)·Vs/dx`, `dx` = element MAX
edge, `Vs = sqrt(mu/rho)` NEAREST-GRID at the element BARYCENTRE.

| mesh | requirement | gate on Vs/dx |
|---|---|---|
| intermediate | 0.5 Hz at p3 | 0.6667 |
| heavy | 0.5 Hz at p3 | 0.6667 |
| heavy | **1 Hz at p5** | **0.8000** ← binding |

## 1. What the gate is judged against — the finding that reshaped this work

The gate had been read from the deck's `safs_material_cvm.nc`: MUSCAL resampled
onto 1500 m lateral / **250 m uniform vertical**. MUSCAL is 0.01 deg laterally
with a **50 m depth step through the top 500 m**. Both agree at the surface;
the deck file simply has no level between 0 and 250 m.

| depth | deck nc Vs p50 | MUSCAL Vs p50 |
|---|---|---|
| 0 | 526 | 534 |
| 50 | *(no level)* | 996 |
| 100 | *(no level)* | 1,429 |
| 200 | *(no level)* | 2,091 |
| 250 | 1,950 | 2,353 |

Nearest-grid on the deck file therefore hands **every barycentre in the top
125 m the z = 0 value**, producing a median **3.78x** step (p90 4.93x, max
8.99x) where MUSCAL has a smooth ramp.

Re-scoring the shipped heavy mesh's 366,083 deck-nc failures against MUSCAL
(`rescore_muscal.py`):

* **only 10.4 % still fail at 0.6667, and 28.1 % at 0.8**;
* 69.8 % of them sit inside the deck file's z = 0 bin, where MUSCAL reports
  **2.18x faster** rock (median 345 -> 751 m/s).

About three quarters of the apparent failures were a **vertical binning
artefact**. The cost difference is not marginal — a fully compliant collar
prices at (`field3d.py`, gradient-limited, x3.66 measured build factor):

| velocity source | hgrad 1.3 | 1.5 | 2.0 |
|---|---|---|---|
| deck nc, gate 0.8 | 28.6M | 21.7M | 16.9M |
| deck nc, gate 0.6667 | 15.9M | 12.3M | 9.8M |
| **MUSCAL, gate 0.8** | 4.4M | **4.0M** | 3.8M |

against the 13.4M the shipped collar already has. Gating on MUSCAL is not only
more honest, it is **cheaper than what is already deployed**.

The forbidden-zone structure also collapses. Under the deck nc the feasible cell
sizes of a slow column split into two disjoint branches, and gradient-limiting
the feasible ceiling costs **+97 %**; under MUSCAL it costs **+10 %**.

> **This is only valid if the deck runs on native MUSCAL material.** The
> `..._MUSCALNATIVE5NC` decks replace the single 250 m cube with five ASAGI
> blocks on MUSCAL's own depth ladder (dz 50 m in the top 3 km), zero vertical
> resampling. **Pair these meshes with that stack.** Running them against the
> old 250 m cube would leave SeisSol propagating through 526 m/s surface
> material the mesh does not resolve.

## 2. Results

### Intermediate — 0.5 Hz @ p3

`results/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal_final.puml.h5`

| | shipped | **new** |
|---|---:|---:|
| tets | 29,385,401 | **34,740,308** (+5,354,907, **+18.22 %**) |
| verts | 5,458,404 | 6,638,003 |
| gate 0.6667 failures (MUSCAL) | 39,281 (0.134 %) | **246 (0.0007 %)** |
| worst Vs/dx | 0.1275 | **0.3101** |
| inverted tets | 0 | **0** |
| fault triangles | 320,560 | **320,560 identical** |
| fault area | 13073.313420396 km2 | **delta exactly 0.000e+00 m2** |
| free surface | 344,862.000 km2, flat at z=0 | **unchanged, still exactly flat** |

**160x fewer failures.** `check_fault_identity.py`: ALL CHECKS PASS (F1 multiset,
F2 area, F3 BC round-trip, F4 flat lid, F5 no inverted).

### Heavy — 1 Hz @ p5 (binding) and 0.5 Hz @ p3

`results/safalt_fb200_1Hz_p5_shakeoutbox_muscal_final.puml.h5`

Two changes: the **parent was swapped**, then the merged mesh was refined.

| | shipped | **new** |
|---|---:|---:|
| parent | `safalt_fb200_deep40km_refine2` 122,162,105 | `safalt_fb200_deep40km_1Hz_p5` **133,699,789** |
| collar | shared 13,405,498 | rebuilt **13,461,912** |
| tets | 135,567,603 | **157,721,257** (+22,153,654, **+16.34 %**) |
| verts | 23,168,926 | 27,677,926 |
| **gate 0.8 failures (1 Hz @ p5)** | 366,083 on the deck nc | **486 (0.0003 %)** on MUSCAL |
| gate 0.6667 failures (0.5 Hz @ p3) | 42,228 on the deck nc | **246 (0.0002 %)** on MUSCAL |
| worst Vs/dx | 0.0868 | **0.3623** |
| inverted tets | 0 | **0** |
| fault triangles | 5,128,960 | **5,128,960 identical** |
| fault area | 13073.313420396 km2 | **delta exactly 0.000e+00 m2** |
| free surface | flat at z=0 | **unchanged, still exactly flat** |

After the weld the merged mesh measured 93,213 failures at 0.8 (parent 7,648 +
collar 85,565); refinement took that to **486**, a 192x reduction.
`check_fault_identity.py`: ALL CHECKS PASS.

| heavy | hops | rounds | tets | failures @ 0.8 |
|---|---:|---:|---:|---:|
| merged | — | — | 147,161,701 | 93,213 |
| pass 1 | 5 | 55 | 157,332,330 | ~7,462 |
| pass 2 | 12 | 80 | **157,721,257** | **486** |

**The shipped heavy was welded to the wrong parent.** `refine2` was only ever
certified at 0.5 Hz/p3 and fails the 1 Hz gate on **274,299 of its own cells** —
75 % of the product's total. Its 1 Hz sibling already existed with those closed
(0 failures at both gates on the deck nc).

It could not reuse the shipped collar. Measured (`wall_compare.py`):

| parent | vertical wall tris | wall verts | canonical sha1 |
|---|---:|---:|---|
| `safalt_0d5Hz_p3_deep40km` (intermediate) | 26,704 | 15,205 | `90c575647b1e5452` |
| `safalt_fb200_deep40km_refine2` | 26,704 | 15,205 | `90c575647b1e5452` |
| `safalt_fb200_deep40km_1Hz_p5` | **46,041** | **25,069** | `24976a75ca1383b0` |

The 1 Hz parent's LEB pass bisected the wall, so a new collar was built against
it: 13,461,912 tets, eta_min 0.0689, **eta<0.05 = 0**, min edge 50 m, wall
vertices reproduced to **1.82e-12 m**. Weld: **46,041 interface faces identical
on both sides**, hull identity OK, 0 negative-volume cells.

## 3. Method, and the knob that mattered

Conforming 3-D longest-edge (Rivara) bisection on the MERGED mesh, fault edges
frozen. Bisection only splits, so compliance is monotone, geometry is exact (a
flat lid stays exactly flat), and it is deterministic — none of which mmg offers
here, since its ±41 % acceptance band licenses coarsening and has twice driven
this gate BACKWARDS.

**The refinement target is the gate itself.** A pooled lower bound exists to stop
a treadmill, and the treadmill is a property of the deck file's binning: with a
3.78x step at the z = -125 m bin edge, one bisection can drop a child's Vs ~4x,
so a pass mints failures as fast as it fixes them. MUSCAL has no cliff — its
50 m ladder means Vs **bottoms out** at the surface rather than falling off it —
so the sequence converges unaided. Measured on the intermediate mesh:

| target | cells flagged | round time | behaviour |
|---|---:|---:|---|
| `bbox` (strict lower bound) | 1,456,239 | 450–885 s | marked flat, ~+24M tets projected |
| `half` (next generation) | 587,899 | 230 s | marked flat at ~590k |
| **`gate`** | **39,281** | **25 s** | peaked at 115k, fell, **+18.2 % tets** |

against 39,281 cells that actually fail. Both pooled variants force the whole
collar surface to ~800 m for no measured benefit. Their logs are kept
(`close_intermediate_HALFPOOL_STALLED.log`) — per the campaign rule that the
runs which failed are what let you tell a wrong metric from a wrong tie-break.

Multi-pass, reseeding a fresh patch on the residual each time (pass 1 stops when
its k-hop patch RIM starts freezing terminal edges — pass 2 ran `froze rim 0`):

| intermediate | hops | rounds | tets | failures |
|---|---:|---:|---:|---:|
| in | — | — | 29,385,401 | 39,281 |
| pass 1 | 4 | 40 | 34,127,870 | ~24,783 |
| pass 2 | 10 | 60 | 34,650,302 | 886 |
| pass 3 | 14 | 90 | **34,740,308** | **246** |

## 4. Regressions and things that got worse

* **The intermediate mesh does not resolve 1 Hz at p5, and got worse at it**:
  gate 0.8 failures went 121,931 -> 455,365 (1.31 %). This is real and expected,
  not a build defect — refining a cell to satisfy 0.6667 moves its barycentre
  toward the surface, which is the same self-reference described in §1, and the
  stricter gate then catches it. The intermediate is only required at 0.5 Hz/p3.
  Use the heavy mesh for 1 Hz.
* **Neither product keeps its parent block bit-identical.** The frozen-parent
  property of the original collar build is deliberately given up: both parents
  carry their own residual on MUSCAL (intermediate 1,160; 1 Hz p5 parent 7,648
  at 0.8), so refining after the weld is the only way to close them. Fault
  identity, free-surface flatness and BC round-trip ARE preserved and checked.
* **Tet counts rise** — intermediate +18.2 %. Quote the LTS cost, not the tet
  count: the added cells are shallow far-field and join COARSE clusters. **Not
  yet measured here** — do it before quoting node-hours.
* **`eta` median** drifts slightly under bisection (measured -0.9 % on the
  earlier 1 Hz p5 build); `eta_min`, `eta<0.05`, `eta<0.1` and min edge were
  bit-identical there. Not re-measured on these products.

## 5. Residual — what it is, and why it is not a floor

Both residuals are the convergence tail, not a structural class, and both sit in
the same place — the corners of the enlarged box, **beyond both velocity models'
data coverage**. MUSCAL's valid box ends at lon -114.02 and the deck nc's grid at
E 705 km, so Vs there is an edge-clamp value, not a measurement.

| | cells | on free surface | median refine still needed | ideal bill | where |
|---|---:|---:|---:|---:|---|
| intermediate @ 0.6667 | 246 | 99.7 % | 1.10x | ~2,000 cells | E 783–786 km (far SE/E) |
| heavy @ 0.8 | 486 | 100 % | 1.20x | ~1,208 cells | E 770–775 km, N 3,976–3,977 km (far NE) |

Refining further chases an extrapolation. Each additional pass costs roughly 10x
the rounds for ~4x fewer cells.

No fault-edge-pinned or wall-pinned exempt class was needed: pass 2/3 ran with
`froze rim 0` and only 3 fault-frozen terminal edges throughout.

## 6. Verification performed

| check | intermediate | heavy |
|---|---|---|
| gate census, MUSCAL, both gates | PASS 246 @ 0.6667 | PASS **486 @ 0.8**, 246 @ 0.6667 |
| fault triangle multiset | PASS 320,560 | PASS **5,128,960** |
| fault area delta == 0 | PASS 0.000e+00 m2 | PASS **0.000e+00 m2** |
| BC round-trip (interior x2, boundary x1) | PASS | PASS 2,564,480 faces |
| free surface flat + area | PASS 344,862.000 km2 | PASS, z in [-0, 0] |
| inverted tets | PASS 0 | PASS 0 |
| weld interface identical | n/a | PASS 46,041 faces |

**Not yet done** (next steps, in priority order): LTS cost / `dt_min` and
`r_insphere` comparison against the shipped meshes — min EDGE is the wrong dt
test, a flat tet keeps an ordinary edge while its insphere collapses; Stage F
deck compatibility (stress/friction nc coverage, hypocentre snap, receiver
containment under the LOCAL free-surface triangle); quality comparison
(`compare_quality.py`) against the shipped products.
