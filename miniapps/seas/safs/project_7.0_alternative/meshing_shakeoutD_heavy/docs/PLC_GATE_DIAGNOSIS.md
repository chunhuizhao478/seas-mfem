# ALT ShakeOut-D: `recoversubfaces` failure -- cause and fix

The assembled ALT PLC (4,603,808 facets) failed boundary recovery with

```
Recovering boundaries...
RuntimeError: Internal TetGen error within `recoversubfaces`.
```

under both `tetgen -d` and the real `-pY`. **Cause: tetgen's boundary recovery is
order-dependent, and this PLC's input ordering is unrecoverable. Fix: permute the
vertex/facet order (`fill_domain.py --shuffle 1`).** Same points, same facets, same
geometry -- so the fix costs nothing.

Resolved. Stage 2 output:

```
9,613,853 tets, 2,300,008 verts        (55 s)
PLC vertices preserved to 2.911e-11 m
volume 12,480,379 km3                  (domain 12,480,379 km3 -- exact)
fault facets surviving as tet faces: 2,564,474 of 2,564,480 -> 6 lost (0.00023 %)
```

## What the `-d` gate does and does not tell you

`isolate_selfint.py` runs bare `switches="d"`. `-d` only *detects*; it never emits a
tetrahedralization, so the wrapper always raises. PREFERRED's final known-good PLC
logged the same pair of lines (`selfint7.log`), so

> "The input surface mesh is correct." + "Failed to tetrahedralize" **is the PASS
> signature.**

ALT's `fault only` and `hull only` were therefore clean. Only `hull + fault` was
different, and it crashed *before* reaching a verdict.

**A correction worth recording:** an earlier revision of this document concluded the
crash was a `-d`-only robustness artifact and that the real fill would succeed. That
was wrong -- `-pY` failed identically. The error was reasoning from local evidence:
every reduced `-pY` box that passed had its lid at z = -0.5, so all of them excluded
the trace, which is the only place the fault and hull meet. Passing boxes that omit
the interaction under test prove nothing about it.

## Structural hypotheses -- all tested, all refuted

None of these was the cause. Recorded so they are not re-litigated.

| hypothesis | test | result |
|---|---|---|
| bad fault/hull weld, near-duplicate vertices | KD weld audit | 8,137 welded at **0.000e+00 m**, all distinct; **0** pairs < 0.5 m in the PLC |
| non-manifold (branch) edges | edge valence census | ALT **840**, PREFERRED **1,048** and it filled |
| T-junctions at the surface | junction verts at z=0 | ALT **2**, PREFERRED **3** and it filled |
| T-junction (vertex on a non-incident edge) | exhaustive local scan | **0** within 1 m |
| needle/sliver triangles (min edge 2.34 m) | reduced `-pY` boxes | recovered fine, 9-52 interior Steiner points |
| shallow "flat wedge" at the trace | tris with >=2 verts at z=0, abs(nz)>0.5 | ALT **232** (max 0.72), PREFERRED **2,340** (max 0.98) |
| fault poking above the free surface | z audit of unwelded fault verts | max **-5.076 m** |
| trace edges never embedded (PREFERRED's real defect) | fault edge set vs top-surface edge set | **0 of 8,137 missing**; 0 *interior* in-plane edges (PREFERRED had 4) |
| sheet-to-sheet proximity at the SE tip | vertex-to-non-incident-triangle scan | **0** pairs under 5 m |
| top-surface quality at the SE tip | triangle quality census | min q **0.644**, min edge 62 m |

## Bisection: where recovery actually breaks

`fill_bisect.py` keeps the hull whole and admits only a subset of the fault.

| admitted fault | facets | result |
|---|---|---|
| everything below -1000 m | 2,391,112 | PASS 9,634,801 tets |
| everything below -100 m | 2,522,535 | PASS |
| everything below -1 m (drops the 15,450 surface-touching tris) | 2,549,030 | PASS |
| + surface band s in [0, 260] km | +4,740 | PASS |
| + surface band s in [260, 455] km | +10,643 | PASS |
| + surface band s in [455, 520] km only | +67 | PASS |
| + surface band s in [260, 520] km | +10,710 | **FAIL** |
| full band minus the 67 tip tris | 2,564,413 | PASS 9,639,936 tets |
| full band minus 88 tris elsewhere (**control**) | 2,564,392 | **FAIL** |

The control is what makes this a localisation rather than a coincidence: removing a
comparable set away from the tip does *not* help.

Those 67 triangles sit at the fault's SE lateral tip (E 616,573-618,343,
N 3,691,705-3,693,687, z -146.27..0), where trace chain `c3` terminates in
mid-surface. They are **not** degenerate: min edge 80.2 m, min quality 0.768,
nearly vertical (abs(nz) 0.02-0.07), 34 trace vertices at ~80 m spacing. Nothing is
wrong with them -- they are simply where the insertion history runs out of flips.

## Why the permutation, and not the alternatives

| approach | outcome |
|---|---|
| input order, `-pY` | FAIL |
| **permuted order, `-pY` (seed 1)** | **PASS -- 279 Steiner, 6 fault facets lost** |
| permuted order, `-pY` (seed 2) | PASS -- 297 Steiner, 8 lost |
| drop `-Y` | PASS -- 1 Steiner, but **215,780 fault facets lost (8.4 %)** |

Dropping `-Y` looks attractive (it adds almost no points) and is the wrong choice:
without it tetgen may re-diagonalise input facets, silently retriangulating 8.4 % of
the frozen fault. The surface is unchanged but the DR facet layout is not, which is
exactly what "reuse the fault triangulation verbatim" forbids.

For scale, PREFERRED's own base fill added **969** Steiner points and lost **2**
facets, so ALT's 279 / 6 is the same order of residue, not a new compromise.

## Note on ZTOP = -60

The deployed ALT mesh has 7,785 fault-boundary vertices that already daylight plus a
390-vertex buried shallow tip (-5.03 .. -142.46 m); exactly **352** of those sit
above -60, which is precisely what `extract_fault_surface.py` snapped
(8,137 - 7,785 = 352). ZTOP=-60 therefore daylights part of the tip the original
build left blind. It is not a build blocker -- the resulting wedges are 10x fewer
and less extreme than PREFERRED's -- and it is not related to the recovery failure,
which sits 60+ km away at the SE terminus. ZTOP=-60 stands as chosen.
