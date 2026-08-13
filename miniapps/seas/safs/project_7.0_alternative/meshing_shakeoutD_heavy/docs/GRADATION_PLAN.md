# ALT: reaching g = 0.15 without a second mmg pass

## The bug

The user specified gradation **g = 0.15**. mmg's `-hgrad` is the RATIO `1+g`, so
that is `-hgrad 1.15`. Both the ALT and PREFERRED builds ran `-hgrad 1.3`, i.e.
**g = 0.30 -- half as gradual as specified**.

## Why this does NOT require re-running mmg

Gradation is a property of the size **field**, not of the mesher. A field that
already satisfies `h(y) <= h(x) + g|x-y|` everywhere is g-graded no matter what
produced the mesh being refined against it. So the fix is to pre-grade the field
at g = 0.15 on a background grid and drive **LEB refinement** of the existing
mmg output against it, rather than to redo a 9 h mmg pass at `-hgrad 1.15`.

This also merges two stages that were separate: the gate term is already in the
field, so **gate closure and gradation happen in the same refinement pass**.

Consequence for ALT: the stage-1 mmg running at `-hgrad 1.3` is a **valid base**.
It is not wasted and does not need repeating.

Method: a min-plus (distance-transform) grading of the field on a background grid.

### The distance metric is NOT a free choice -- L1 under-refines

A separable forward+backward sweep per axis with additive cost computes the **L1**
transform, and an earlier revision of this document claimed L1 was "conservative --
it can only make the field smaller, never larger, so it cannot under-refine."
**That is backwards.** The transform is

```
h(x) = min over y of [ h0(y) + g * d(x,y) ]
```

and since `L1 >= L2`, every candidate `h0(y) + g*d_L1` is >= its Euclidean
counterpart, so the minimum is too:

```
h_L1(x)  >=  h_euclid(x)      pointwise -- L1 grading is COARSER
```

Equivalently `h(x) <= h(y) + g*d_L1` is a *weaker* constraint than the Euclidean
one, because the bound is larger. Measured, single 115 m seed, g = 0.15:

| offset | L2 | L1 | h_euclid | h_L1 | L1 coarser by |
|---|---|---|---|---|---|
| (100,0,0) | 100.0 | 100.0 | 130.0 | 130.0 | 0.0 % |
| (100,100,0) | 141.4 | 200.0 | 136.2 | 145.0 | 6.5 % |
| (100,100,100) | 173.2 | 300.0 | 141.0 | 160.0 | 13.5 % |
| (500,500,500) | 866.0 | 1500.0 | 244.9 | 340.0 | 38.8 % |
| (2000,2000,2000) | 3464.1 | 6000.0 | 634.6 | 1015.0 | 59.9 % |

A field that is g-Lipschitz in L1 is only `g*(|u|_1/|u|_2)`-Lipschitz in Euclidean,
so an L1 sweep at g = 0.15 actually delivers

| direction | \|u\|_1/\|u\|_2 | g_eff |
|---|---|---|
| axis | 1.000 | 0.150 |
| face diagonal | 1.414 | 0.212 |
| body diagonal | 1.732 | **0.260** |

Exact on the grid axes, worst on the body diagonal -- the generic direction away
from a dipping fault. Since the entire point of this change is to deliver the
specified g = 0.15 rather than 0.30, an L1 sweep would land at 0.15-0.26 and, worse,
do so **anisotropically**: transitions look gradual along the grid axes and steep
along diagonals, which is a nastier artefact than a uniform 0.30.

Options, cheapest first:

1. **Sweep at `g/sqrt(3)` = 0.0866.** Then `g_eff <= 0.15` in every direction. One
   line, keeps the exact separable transform, but over-refines on the axes by up to
   sqrt(3) and so costs more than the 110.5 M predicted for true Euclidean g=0.15.
2. **3-D chamfer (Borgefors) mask** -- 6 face neighbours at 1, 12 edge at sqrt(2),
   8 corner at sqrt(3), each weighted by the REAL anisotropic spacings.
   Approximates Euclidean to a few percent instead of 73 %. **Preferred.**
3. Felzenszwalb-Huttenlocher is separable and exact for *squared* Euclidean, but our
   cost is linear `g*d`, so it does not drop straight in.

### Built and VERIFIED for ALT

```
grid    1192 x 1138 x 161 = 218,395,856 nodes (isotropic 500 m)
fault   2,564,480 facets -> 256,448 decimated centroids
grade   Euclidean min-plus, R=2 (98 offsets), sweep g = 0.14286 (= 0.15/1.05)
        converged in 66 iterations (sum drop -> 0)
h       115 .. 5000 m          out: build_tmp/target_g015.npz, 0.81 GB
```

The point of the whole exercise is the *delivered* gradation, so it is measured
rather than assumed -- and measured **on the body diagonal**, which is where an
L1 or anisotropic-grid field fails while still looking correct on the axes:

| direction | max abs(dh)/d | 99.99 pct |
|---|---|---|
| x | 0.14286 | 0.14286 |
| y | 0.14286 | 0.14286 |
| z | 0.14286 | 0.14286 |
| **body diagonal** | **0.14286** | 0.14286 |

Uniform at the swept value in every direction, so the delivered Euclidean
gradation is <= 0.143 < 0.15 everywhere. For contrast, an L1 sweep at g=0.15 on
the natural non-uniform z would have read ~0.26 on that last row while showing a
correct 0.15 on the first three -- the failure mode is invisible unless the
diagonal is checked.

### The z axis is non-uniform

The background grid uses 100 m spacing near the surface, 500 m at mid-depth and
2000 m deep. The per-axis sweep must use the **actual level spacing**, not a scalar
`dz`: a constant would over-grade the deep half and under-grade the top 3 km, which
is exactly where the gate bites. This compounds with the metric choice above -- a
chamfer mask needs the true anisotropic spacing per neighbour offset.

Because the field is g-graded by construction, the pooled lower bound needed by
the refinement loop is **exact** and needs no second lookup:

```
h(y) >= h(x) - g|x-y|      =>      pooled(x, dx, frac) = h(x) - g*frac*dx
```

This also sidesteps the vertical-extent pooling runaway recorded in the project
notes (+505 k tets for a flat failure count).

No change is needed to the LEPP machinery: the loop only asks `at(p)` and
`pooled(v,dx,frac)`, both compared as `Vs/gate`, so a field target returning
`h_target*gate` makes the identical code refine to `dx <= h_target`.

## What g = 0.15 actually costs

An earlier estimate of **140-210 M tets / 29-43 GB ("may need a cluster")** was
wrong and is withdrawn. It applied a "7.13x gradation multiplier" =
70.6 M actual / 9.90 M predicted. That baseline is invalid: the 9.90 M is an
integral of the **floored** metric, and the floored metric is unachievable by
construction -- the fault triangles are frozen at 115 m, so mmg cannot place
1400 m tets against them whatever the `.sol` says. 7.13x measured that
prediction's *error*, not a gradation effect. The full-field integrals already
contain h = 115 m at the fault and grade outward, so multiplying them by 7.13
double-counts the same band.

Falsification: PREFERRED's full-field g=0.30 prediction was 106.0 M and its
actual floored stage-1 produced 70.6 M -- **below** the prediction, as a floor
should be. A real 7.13x inflation cannot put 70.6 M under 106 M.

Measured (size_budget, hmax 5000, gate 0.8, k=18.8992), both trees:

| | g=0.15 | g=0.30 | ratio |
|---|---|---|---|
| ALT full field | 110.5 M | 89.7 M | **1.23x** |
| PREFERRED full field | 123.8 M | 103.0 M | **1.20x** |
| ALT gate only | 43.4 M | 23.6 M | |
| ALT fault treatment | 67.1 M (61 %) | 66.1 M (74 %) | |

Two independent trees agree on ~1.2x, so g = 0.15 costs roughly **17-18 GB** at
the measured 14.5 GB / 70.6 M scaling -- inside 36 GB. No cluster needed.

**The gate is not nearly free.** A separate claim that the 1 Hz-at-p5 gate costs
only 14.64 M, making >93 % of the mesh fault treatment, used the *ungraded*
integral and compared it against a *graded* total. Like for like, ALT's gate
alone is 23.6 M at g=0.30 and 43.4 M at g=0.15, so the fault treatment is
**61-74 %** -- still dominant, but the frequency requirement is ~3x more
expensive than that framing suggests. This matters if the fault-zone treatment is
ever put up for negotiation, which is the one thing the build was told to preserve.

These are k-calibrated integrals, not mmg runs. The single hard datum -- 70.6 M
actual against a 106.0 M full-field prediction -- says they are the right order
and conservative.

## Note on the stage-1 metric floor

`build_metric.py --hmin 1400` floors the metric at **100.0 % of free-surface
vertices and 100.0 % of fault vertices** (95.6 % overall), so mmg is told to
coarsen the gate-correct free surface `surface_field.py` built at 115-1143 m.
PREFERRED measured the effect: top surface 1,829,708 -> 250,839 triangles.

Deliberately **not** restarting for this. The trace rides on the frozen fault
RequiredTriangles so it survives; the top is flat, so re-refinement is
geometrically exact; and PREFERRED's post-stage-1 gate was still only 0.749 %
sub-gate. Stage 2 should floor **only the fault term** and leave the gate term
exact.
