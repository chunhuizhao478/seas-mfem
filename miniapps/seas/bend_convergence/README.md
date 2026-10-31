# Bend convergence test

Chunhui Zhao, 2026-10-30: "prepare the bend convergence test".

- **Question:** at a bend, is the drop of sigma_n' per metre of slip set by
  the elasticity or by the mesh size h?
  - If it grows as h shrinks, it is set by the mesh: refining makes it worse.
  - If it settles, it is a physical limit that a bend treatment must handle.
- **Why it matters:** the UCERF runs reach sigma_n' of about -270 MPa at a
  few bend nodes (the 1 MPa floor only hides it from friction). The earlier
  500 m bend runs gave about 5 MPa/m per degree at a sharp kink, close to
  4 mu theta / h, which would double with each halving of h.

## Runs

Nine runs: three shapes, each at three mesh sizes.

| Shape | Tag | Bend | Earlier 500 m run |
|---|---|---|---|
| sharp | `sharp05` | 5 deg at (0, 0) | stopped at sigma_n' < 0, 3.69 s into the event |
| 5 km arc | `arc05L05` | 5 deg over s = -2.5 to 2.5 km, R = 57.3 km | stopped at 11.61 s (500 m facets, 0.5 deg kinks) |
| 10 km arc | `arc05L10` | 5 deg over s = -5 to 5 km, R = 114.6 km | hand-back passed (0.25 deg kinks) |

- **Mesh sizes:** h = 500, 250 and 125 m on the fault (`_h500`, `_h250`,
  `_h125`).
- **Fault points on the true curve at every h:** the facet angle
  theta h / L_arc halves with h. The earlier arc meshes kept 500 m facets.
- **Model:** plane strain, no sigma_n' floor and no 2.5D term, so the bend
  response is measured on its own and a run stops where sigma_n' reaches 0.
- **Decks:** `config/fd_cycle/mapview_rs_fd.toml` (adf231e6) with the
  earlier bend runs' edits.
  - Physics: rate and state, a = 0.004 for |x| <= 10 km (0.04 from 12 km),
    Dc 0.14 m, sigma_n' 25 MPa.
  - Geometry: a 30 km fault in a +-120 km box, with absorbing walls on every
    side.
  - Run length: at most 2 dynamic phases.
  - Edits:
    - `[mesh] path`;
    - `max_dynamic_phases = 2`;
    - the prefix and `handback_table = true`;
    - `region_x3 = [-1e4, 1e4]`;
    - new here: fault snapshots every 0.01 s in the dynamic regime (was
      0.1 s), so the fit below has points even where sigma_n' reaches 0
      within 0.3 m of slip.
  - The nine decks differ only in the mesh path and the prefix.

## What each outcome means (estimates, not measured)

- **Sharp bend:** a kink in a slipping fault acts like an edge dislocation
  of Burgers vector about theta D at the vertex. The stress there is
  singular, so the vertex value should grow about as 1 / h: about 22, 45
  and 90 MPa/m at 500, 250 and 125 m (from 4 mu theta / h).
- **Arc:** each facet kink of theta h / L_arc gives about
  4 mu (theta h / L_arc) / h = 4 mu theta / L_arc, the same at every h:
  - about 2.2 MPa/m for the 5 km arc;
  - about 1.1 MPa/m for the 10 km arc.
  This is the curvature effect, of order mu D kappa, and it should converge.
- **Reading the result:** if the arcs converge and the sharp bend does not,
  the UCERF bend tension is real for its curvature (R = 10 km arcs: about
  4 mu / R = 13 MPa/m), not a mesh artefact.

## Files

- `mesh/bend.geo`: one parametric mesh (theta_deg, arc, lc_fault).
  - It is `mapview_rs/mesh/mapview_rs.geo` with a +-120 km box and the
    fault built on the curve.
  - Its header has the gmsh commands.
  - Two more inputs for `simple_cases/` (2026-10-31): `box` (the box
    half-width, default 120 km) and `s_bend` (the bend's position along
    the fault, default 0). With the defaults the nine meshes below come
    out byte for byte the same.
- `mesh/*.msh`: the nine meshes, committed (cloud gmsh 4.15.2).
  - They were checked to be on the curve: 60, 120 and 240 fault edges of
    exactly h, facet kinks theta h / L_arc, arc radius error below 1e-10 m.
  - sha256 (also in the run script):

    | Mesh | Triangles | sha256 |
    |---|---|---|
    | sharp05_h500 | 14,760 | be86e9ce5b9740d2dea6c66c8c63527cae59eff70c70e362eedcf18043695c26 |
    | sharp05_h250 | 23,524 | 7610004069c69ef7e89a2e565e40304a4e0e762e8fd3ae83482670dd085de41d |
    | sharp05_h125 | 38,364 | 13b3d0d02930d2e295b1d01791e172e52357759dd04379919313ffd083de1e8c |
    | arc05L05_h500 | 14,752 | 0f8a78206be07cef7ab74f8ca28d823eefc8262a5d3a381ed13b43c13870259d |
    | arc05L05_h250 | 23,484 | c953ea83c280076ff90d0a8181628935995dd0fb1b766d2b6bb37c827dfef795 |
    | arc05L05_h125 | 38,282 | 9ba4d920668bff7691bf76052a9704c1c00ae7017b01b9d338ad62a53163dddc |
    | arc05L10_h500 | 14,854 | e00f734b9159b23f888183274c71e63fd220899c8728aed1bd2b26c3eeabf536 |
    | arc05L10_h250 | 23,586 | 7e9b649e2048a51dd9cacbf5045a77d300cd2d438ad47abf6aedf46a1279f781 |
    | arc05L10_h125 | 38,246 | 42673988d0974546e06da6fc6a433a9ab6924084a2947ad080e1ec1ebe5e2255 |

- `configs/<tag>.toml`: the nine decks.
- `scripts/run_bends_mac.sh`: the nine runs on the Mac, one after the other.
  It checks:
  - the GitHub-first guard;
  - that the binary's code equals HEAD (`BUILD_COMMIT_MAC_FDCYCLE`);
  - the mesh sha256s.
  It takes the Mac run lock and sets OMP_NUM_THREADS=1. A run that stops at
  sigma_n' < 0 is a result, so the script moves on to the next deck.
- `scripts/analyse_bends.py`: the comparison (`summary.txt`, `summary.json`,
  `bend_convergence.pdf`). Per fault point (once per face, so the two sides
  of a kink stay apart), it computes:
  - k, the slope of sigma_n' against slip over the snapshots with slip >=
    0.02 m;
  - per run: the lowest sigma_n' (where, the slip there, the fitted and
    secant k there) and the most negative k along the fault;
  - k at the bend points (vertex; arc start, middle, end);
  - the ratio per halving of h.

## Workflow (Mac, from `miniapps/seas` of the checkout at the pushed commit)

```bash
conda activate mumps-petsc
make -f tests/2d_inplane/2d_inplane.mk seas_fd_cycle_driver
git rev-parse HEAD > BUILD_COMMIT_MAC_FDCYCLE
R=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/bend_convergence
bend_convergence/scripts/run_bends_mac.sh $R 8
python3 bend_convergence/scripts/analyse_bends.py $R
```

- **Run time:** about 1 to 2 h for the nine on 8 ranks (scaled from the
  500 m runs, not measured).
- **Disk:** a few GB (inferred).
- **First check:** `sharp05_h500` should stop near 3.69 s into the event, as
  before (same physics; only the reporting region and the output interval
  changed).

## Checks done before the runs (cloud, 2026-10-30)

- **Meshes:** the nine built, and their geometry checked as above.
- **Decks:** each was diffed against `mapview_rs_fd.toml`. Only the edits
  listed above differ.
- **`run_bends_mac.sh` with a stub mpirun:** it refused:
  - a missing stamp;
  - an unknown tag;
  - an existing output folder (skipped).
  With a stamp it wrote `RUN_INFO.txt` with the commit, mesh sha256 and
  exit code.
- **`analyse_bends.py` on synthetic snapshots in the writer's layout:** the
  inputs were per-face points, with k = -4 mu theta / h at the vertex and
  -4 mu theta / L_arc on the arcs, stopping at sigma_n' < 0.
  - Sharp: k by h came back as -22.4, -44.7 and -89.5 MPa/m (x2.00 per
    halving).
  - Arcs: -2.24 and -1.12 MPa/m at every h (x1.00).
