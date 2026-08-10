# QUALITY REPORT — ALT intermediate / heavy extended to the ShakeOut box

Scope: **ALT intermediate and ALT heavy only.** Parent frozen; only the far
field outside the parent's wall is built. Domain box and its derivation:
`DOMAIN_SHAKEOUTBOX.md` (714.000 x 483.000 km, contains the ShakeOut grid box
with 2.6–16.4 km margin and contains the PREFERRED box, so ALT n PREF is
exactly the PREFERRED box).

## Requirements

| tier | requirement | gate on Vs/dx |
|---|---|---|
| intermediate | 0.5 Hz at p3 | 0.6667 |
| heavy | 0.5 Hz at p3 | 0.6667 |
| heavy | **1 Hz at p5** | **0.8000** ← binding |

Gate convention (locked): f = (p/4)·Vs/dx, dx = element MAX edge, Vs = sqrt(mu/rho)
NEAREST-GRID at the element BARYCENTRE, from the deck's `safs_material_cvm.nc`.

## Current state (uniform 2,500 m far field)

| | intermediate | heavy |
|---|---|---|
| tets | 29,385,401 | 135,567,603 |
| parent | 15,979,903, **bit-identical** | 122,162,105, bit-identical |
| collar | 13,405,498 | same collar |
| weld | 26,704 interface faces identical | identical |
| hull identity | tagged 2,452,624 == hull 2,132,064 + fault 320,560 | OK |
| collar eta_min | 0.0779 | 0.0779 |
| collar eta<0.05 | **0** | 0 |
| collar min edge | 45 m | 45 m |
| gate 0.6667 | 42,228 collar cells fail (0.315 %) | same |
| gate 0.8000 (1 Hz p5) | 91,784 collar cells fail (0.685 %) | same |

## THE MAIN FINDING — the gate is self-referential, and fine is worse

The gate samples Vs at the element BARYCENTRE. A cell of vertical extent h
resting on the free surface has its barycentre at roughly z = -h/3. The CVM's
z = 0 bin is its slowest (Vs 155–545 m/s over this footprint) and the next bin
down is ~4x faster. So **making surface cells smaller drags their barycentres
INTO the slowest bin and tightens the constraint they are being refined to
satisfy.** The refinement manufactures its own demand.

Measured, same footprint, same parent, same tool:

| far-field design | collar tets | gate 0.6667 failures | eta_min |
|---|---|---|---|
| gate-driven + interior seed | 22,266,601 | 1,226,028 (5.506 %) | 0.0663 |
| **uniform 2,500 m** | **13,405,498** | **42,228 (0.315 %)** | **0.0779** |

40 % fewer cells, 29x fewer failures, better quality. The earlier design was
worse on every axis.

Corroborating sweep — uniform h vs volume fraction meeting 0.5 Hz:

| h | cells (analytic) | volume at 0.5 Hz |
|---|---|---|
| 500 m | 362,975,040 | 99.91 % |
| 1,000 m | 45,371,880 | 99.38 % |
| 2,500 m | 2,903,800 | 98.85 % |
| 3,000 m | 1,680,440 | 98.32 % |

500 m costs **216x** the cells of 3,000 m and buys **1.6 percentage points**.
The shortfall is a thin skin at z = 0, not a volume problem.

## NEXT STEP — depth grading, for minimum mesh gain

Uniform 2,500 m is still wrong in the other direction: it over-refines at
depth, where Vs 2,800–3,300 m/s permits 3.5–5 km cells. Solving the
self-consistency h <= Vs(-h/3)/gate on a lateral p10 of Vs gives

| z | Vs(p10) at barycentre | h allowed |
|---|---|---|
| 0 | 2,016 | 3,023 m |
| -3,023 | 2,800 | 4,200 m |
| -12,051 | 3,383 | 5,075 m |
| -28,141 | 4,194 | 6,000 m |

`code/depth_profile.py --gate {0.6667,0.8} --pct N` produces this.

**CAVEAT, measured:** the analytic cell estimate `6·A·slab/h^3` UNDERPREDICTS
the real build by ~4.6x — it predicted 2.9M for uniform 2,500 m and tetgen
produced 13,405,498. Treat the profile's ~420k as ~2M in practice, and always
re-measure after the build rather than trusting the estimate.

**The percentile is the real knob.** p10 sizes the profile to the bulk far
field and leaves the slowest ~10 % of columns short at the surface; the
absolute minimum column (Vs 155 m/s at z = 0, 795 m/s at z = -1 km) cannot be
resolved at any affordable cost and is a geometry/material-forced floor, not a
fixable defect.

## Open items

1. **gate 0.8000 (heavy, 1 Hz p5) MEASURED on the uniform collar**: 91,784 of
   13,405,498 fail (0.685 %), worst 0.0868.  Failures are NOT a surface skin
   here -- median barycentre -920 m, p10 -4,314 m, and 10,959 of them below
   -4 km -- so they are the genuinely slow CVM columns, and 2,500 m is simply
   too coarse for them at this gate.  A heavy-specific collar at gate 0.8
   (h ~ 2,000 m by the same profile) is the remaining build.
2. The 42,228 residual at 0.6667 is now small enough for exact closure by local
   Rivara LEB (skill architecture 4) — it was not, at 1.23M.
3. LTS cost of the uniform collar not yet measured. The gate-driven collar
   measured 2.144x for 2.393x the cells with dt_min unchanged; the uniform
   collar should be cheaper still, but **measure before quoting node-hours**.

## Regressions to note

* the shipped uniform collar has collar min edge 45 m against the earlier
  61 m — still far above the intermediate parent's own 9.34 m, so no new dt
  floor, but it is a regression against the previous collar;
* worst Vs/dx over the collar is 0.0868 (uniform) vs 0.2615 (gate-driven):
  FEWER cells fail, but the single worst cell is worse. Both are the slowest
  CVM columns; neither design resolves them.
