# From-scratch ShakeOut-box HEAVY mesh — design, budget, and build plan

Goal: rebuild the heavy mesh on the ShakeOut domain **without** inheriting the
small-domain configuration, keeping the fault-normal refinement, with gradual
size transitions everywhere, holding **0.5 Hz @ p3 (Vs/dx >= 0.6667)** and
**1 Hz @ p5 (Vs/dx >= 0.8000)**.

Domain: axis-aligned **714.000 x 483.000 km**, z 0 .. -40,000 m, flat top at
z = 0 (`../DOMAIN_SHAKEOUTBOX.md`).

## 1. The spec, MEASURED off the deployed mesh (`code/measure_envelope.py`)

Not re-derived from a remembered "200 m within 1 km" — the deployed heavy is a
base fb200 build plus two red refinements plus a gate campaign, so its real
envelope is not any single documented rule.

| distance from fault | median max edge | p95 |
|---|---:|---:|
| 0–100 m | **93 m** | 129 |
| 100–200 m | 100 | 138 |
| 200–400 m | 112 | 156 |
| 400–700 m | 177 | 227 |
| 700–1,000 m | 218 | 293 |
| 1–1.5 km | 244 | 352 |
| 1.5–2 km | 357 | 530 |
| 2–3 km | 518 | 795 |
| 3–5 km | 761 | 1,280 |
| > 12 km | 759 | 2,690 |

The near-fault band is **~110 m at EVERY depth** (110/110/106/115/107 m from the
surface down to -20 km), i.e. a uniform fault-normal shell, not a depth-tapered
one. Far field is gate-driven: 462 m median above -1 km, ~2,300 m below -6 km.

## 2. The budget — where the 160.8M cells actually are

| band | cells | share |
|---|---:|---:|
| within 1 km of the fault | **95,384,826** | **59.3 %** |
| 1–12 km transition | 12,783,480 | 7.9 % |
| beyond 12 km (pure gate) | 52,657,956 | 32.7 % |
| **total** | **160,826,262** | (deployed: 160,825,939 — 0.0002 % match) |

**A from-scratch build to this same envelope is the same ~160M cells.** The
rebuild does not make the mesh smaller; it changes only how it was assembled.

## 3. Is the grading actually bad? MEASURED — no (`code/grading_check.py`)

Adjacent-cell max-edge ratio across every interior face (intermediate mesh,
75,986,963 interior faces):

| | p50 | p90 | p99 | max | > 2.0 |
|---|---:|---:|---:|---:|---:|
| all interior faces | 1.055 | 1.281 | 1.602 | 6.321 | 0.032 % |
| **weld (parent ↔ collar)** | 1.056 | 1.234 | **1.497** | **2.846** | 0.024 % |
| non-weld interior | — | — | 1.605 | **6.321** | — |

**The frozen-parent weld is BETTER graded than the mesh average**, and the worst
transition in the mesh (6.32x) is in the interior, nowhere near the seam. So the
motivation for a rebuild — an abrupt size change at the domain join — is not
present in the deployed mesh. Whatever a rebuild is worth, it is not that.

## 4. Feasibility: a single-shot 160M fill is NOT possible here

This machine has 36 GB. tetgen needs roughly 150–250 B/tet including its
neighbour structures (~25–40 GB at 160M, before the PLC and the point cloud);
gmsh's 3-D Delaunay is several times worse. The 13.4M-cell collar fill took
110 s and a few GB — that scales to a fill that does not fit.

So the build MUST be staged, which is also how this lineage was actually
produced (`safalt_fb200` -> `refine2` = fault edge / 4):

| stage | what | cells |
|---|---|---:|
| A | PLC: ShakeOut box + flat lid + CFM fault strands, conforming curves | — |
| B | surface conditioning: fault at ~440 m, LEB refine + feature-aware collapse | — |
| C | **base volume fill**, band /4 (440 m) and far field /2, graded field | **~9.7M** |
| D | mmg de-sliver, flips, reorient; SAFS tags; PUML `--validate` | ~9.7M |
| E | fault-band refinement back to ~110 m (two red/LEB levels) | ~+94M |
| F | gate closure to 0 at 0.8, `leb_muscal.py --pool zpool --pool-radius 4` | ~160M |

Stage C at ~9.7M is a single comfortable tetgen run. Stages E–F use the tooling
already proven this session on 160M-cell meshes at 10–21 GB.

## 5. Status

Stages 1–4 above (design measurement, budget, grading verdict, feasibility) are
DONE and are what this folder currently contains. The build stages A–F are
specified but NOT yet run.

Tools here:

* `code/measure_envelope.py` — recover the deployed size-vs-fault-distance design
* `code/grading_check.py` — adjacent-cell size ratios, split by weld vs interior
* `code/price_build.py` — design field `min(band envelope, gate ceiling)`,
  Lipschitz-limited, integrated. **Caveat:** its KD distance field over a
  52M-node lattice with a 24 km `distance_upper_bound` did not finish in 35 min;
  the budget in §2 was taken from the direct envelope census instead, which is
  both exact and cheap. Shrink the radius or the lattice before reusing it.

## 6. Recommendation before spending stages A–F

Two measured facts are worth weighing first:

1. the rebuild reproduces **the same ~160M cells** — it is not a cost saving;
2. the weld it would remove is **not a grading defect** (p99 1.497, max 2.846,
   better than the mesh interior).

The deployed mesh already holds both gates exactly (0 failures at 0.8, worst
landing on 0.8000), with the fault triangulation bit-identical to its parent.
A rebuild buys a clean single-lineage provenance and a field-designed grading
everywhere, at the cost of re-running Stage A–F and re-verifying deck
compatibility from scratch. That is a real benefit for a mesh meant to outlive
this campaign — but it should be chosen deliberately, not as a fix for a seam
problem that measurement says does not exist.
