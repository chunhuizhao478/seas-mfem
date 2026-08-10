# The ShakeOut domain box used by the ALT meshes

All three `meshing_shakeoutbox_*` products share ONE box:

    E  72,000 .. 786,000        714.000 km
    N 3,524,000 .. 4,007,000    483.000 km      = 344,862 km2

## How it was derived (same method as the PREFERRED build)

ShakeOut v1 is USGS ShakeMap `sclegacyshakeout2full_se`. Its `grid.xml`
declares lon -121.5 .. -114.0, lat 32.0 .. 36.0 at 0.02 deg (376 x 201 nodes).
Projected to UTM 11N (EPSG:32611, the SAFS CRS) by a **dense sample of the box
boundary** -- 2,001 points per edge, not the four corners, because the
projected edges bulge -- that footprint is

    E  74,758.0 ..  783,423.2    708.665 km
    N 3,540,435.7 .. 3,993,324.3  452.889 km

reproduced here to **0.0 m** against the independent PREFERRED derivation.
Corners alone would understate the SOUTH edge by 3,934 m.

Four requirements then fix the domain box:

1. contain that footprint with about one far-field cell of margin
   (W 2,758 / E 2,577 / S 16,436 / N 13,676 m);
2. contain the parent -- a frozen-parent collar can only ADD;
3. leave a workable collar width everywhere.  The ALT parent's north corner is
   at N 3,996,866.9, i.e. 3,543 m NORTH of the ShakeOut footprint, so the north
   wall is pushed 10,133 m above it.  A wall drawn at the corner itself pinches
   the collar to zero width, and gmsh then reports `2 intersections in the 1D
   mesh` and emits **no elements at all**;
4. **contain the PREFERRED domain box** (E 72,000..786,000,
   N 3,524,000..3,996,000).

Snapped to round km.

## Why requirement 4 matters

With it, **ALT n PREFERRED is exactly the PREFERRED box** (714.000 x 472.000
km), so the common ALT-vs-PREF-vs-ShakeOut comparison mask inherits PREF's
margins on all four sides:

| | W | E | S | N |
|---|---|---|---|---|
| margin of ALT n PREF vs the ShakeOut footprint | 2,758 | 2,577 | 16,436 | 2,676 m |

## What this replaced, and why

The first ALT build used the **contour product's** bbox from
`~/Downloads/ShakeOut2008/04_comparison/HEAVY/shakeout_comparison_HEAVY.json`
(E 74,850.1..781,468.9, N 3,542,641.5..3,993,303.4) -- where the contour data
actually is, which is smaller than the declared grid. Two consequences, both
measured:

* the absorbing wall sat **exactly on the comparison area** with **zero
  margin** on west, east and south.  PGV within about a cell of an absorbing
  boundary is not trustworthy;
* the box did not contain the declared grid box (short 92 m W, 1,954 m E,
  2,206 m S), so it failed that containment test.

Neither old box could simply adopt the other's: the two parents protrude from
each other's boxes on OPPOSITE sides (ALT 10,867 m north of PREF's box; PREF
9,956 m south of ALT's), and a frozen-parent collar can only add. Rebuilding
ALT to the union was therefore the only way to get a common comparison domain
without touching PREFERRED.

## Remaining differences from PREFERRED (not resolvable by a collar)

| | ALT | PREFERRED |
|---|---|---|
| depth | 0 .. -40,000 m (small tier: -22,100.9) | +25.8 .. -39,329.4 m |
| lid | **exactly** flat at z = 0 | flat to +-25.8 m |

Both are inherited from the parents.
