# Simple cases (2D plane strain)

Chunhui Zhao, 2026-10-31: "we will return back from 2.5D to 2D plane
strain"; "I would like to test in a few simple settings"; "go ahead with
the simple cases, run them on local Mac".

- **Question:** the UCERF 2D runs have two problems: tension at fault
  bends, and the dynamic-to-QD hand-back. These runs split them:
  - Does the trouble need a bend? (straight against bent)
  - Does it need the dynamic regime? (QD only against fully dynamic)
- **Plan:** the Claude doc "2D plane strain simple tests plan"
  (https://claude.ai/artifact/TXr2e3UVfsxqUjXKdroUqm): test matrix, mesh
  figure, what each outcome means.

## Runs

Twelve runs, listed in run order (D1f, D2f and Q1 answer the main
questions, so they go first).

| Case | Mesh | Regime | Floor | Why |
|---|---|---|---|---|
| D1f_sharp015_floor | sharp015_h500 | fully dynamic | 1 MPa | D1 with UCERF's floor: does a bend break the hand-back? |
| D2f_sharp05_floor | sharp05_h500 | fully dynamic | 1 MPa | the same at 5 deg |
| Q1_sharp015 | sharp015_h500 | QD only | none | Andy's [2]: UCERF's 1.5 deg angle |
| Q0_straight | straight_h500 | QD only | none | Andy's [1]: the control |
| Q2_sharp05 | sharp05_h500 | QD only | none | Andy's [2]: 5 deg sharp |
| Q3_arc05L10 | arc05L10_h500 | QD only | none | Andy's [2]: 5 deg over a 10 km arc |
| Q4_sharp015vs | sharp015vs_h500 | QD only | none | 1.5 deg in the creeping part (like the UCERF NaN node) |
| D0_straight | straight_h500 | fully dynamic | none | Andy's [3]: the control |
| D1_sharp015 | sharp015_h500 | fully dynamic | none | Andy's [4]: 1.5 deg |
| D2_sharp05 | sharp05_h500 | fully dynamic | none | Andy's [4]: 5 deg sharp |
| D3_arc05L10 | arc05L10_h500 | fully dynamic | none | Andy's [4]: 5 deg arc |
| D4_sharp015vs | sharp015vs_h500 | fully dynamic | none | 1.5 deg in the creeping part |

- **No floor (all but D1f, D2f):** a run that reaches sigma_n' <= 0 stops
  there, and that time and place is a result.
- **Bent runs without a floor (D1, D2):** they will likely stop in event 1,
  before any hand-back (estimate). D1f and D2f test the hand-back on a
  bend instead.

## Setup

All twelve decks are `config/fd_cycle/mapview_rs_fd.toml` with the edits
below. `scripts/make_decks.py` makes them, and the run script checks that
the committed decks are exactly its output.

- **Physics (unchanged):** a 30 km rate-and-state fault; a = 0.004 for
  |x| <= 10 km rising to 0.04 at 12 km; b 0.03, Dc 0.14 m, sigma_n' 25 MPa;
  plate loading at V_p = 1e-9 m/s through the box sides; switch at V_on
  0.01 m/s, back at V_off 0.001 m/s after 10 s; all walls absorbing in the
  dynamic regime.
- **Edits, every case:**
  - `[mesh] path`: `simple_cases/mesh/<mesh>.msh` (+-240 km box).
  - `[time] tfinal`: 12,000 yr (was 1,500 yr). See "Length" below.
  - `[cycle.output] prefix` = the case; the fault stations (below).
  - `[cycle.benchmark] region_x3 = [-1e4, 1e4]`: the bent limb reaches
    y = 1.31 km (as the bend convergence decks).
  - `[cycle.snapshots.fault] dt_dyn = 0.05 s` (was 0.1 s): sigma_n'
    against slip at the bend. Output only.
  - `[cycle.snapshots.volume] dt_dyn = 1.0 s` (was 0.5 s): 1 MB per
    snapshot, so about 0.27 GB per fully dynamic run of 3 events. Output
    only.
- **QD only (Q):** `[cycle] switching = false`.
- **Fully dynamic (D):** `[cycle.output] handback_table = true`.
- **Floor (D1f, D2f):** `[cycle] sigma_n_floor_pa = 1.0e6` (UCERF run8's).

### Stations

All on the fault trace, placed by arc length u from the fault centre
(u = x on the straight west part), so the stations on a bent limb lie on it:

- `fltst_s-12` ... `fltst_s+12`: every 3 km of u.
- `bend_w`, `bend_e`: 1 m from the bend vertex (the arc's mid point) along
  the west and the east face. Each face keeps its own sigma_n' at a kink;
  a station exactly at the vertex would take one face only.
- `bend_s-3.0`, `bend_s-1.0`, `bend_s+1.0`, `bend_s+3.0`: along the fault
  from the bend (Q4/D4: `bend_s+2.0` instead of +3, since the tip is 2.5 km
  past that bend; the arc: also its ends `bend_s-5.0`, `bend_s+5.0`).

## Meshes

Made in the cloud with gmsh 4.15.2 from `../bend_convergence/mesh/bend.geo`
(its new inputs `box` and `s_bend`; its defaults still give the nine bend
meshes byte for byte):

```bash
cd simple_cases/mesh
G=../../bend_convergence/mesh/bend.geo
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 0   -setnumber arc 0     -setnumber s_bend 0     -setnumber lc_fault 500 -o straight_h500.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 1.5 -setnumber arc 0     -setnumber s_bend 0     -setnumber lc_fault 500 -o sharp015_h500.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 5   -setnumber arc 0     -setnumber s_bend 0     -setnumber lc_fault 500 -o sharp05_h500.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 5   -setnumber arc 10000 -setnumber s_bend 0     -setnumber lc_fault 500 -o arc05L10_h500.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 1.5 -setnumber arc 0     -setnumber s_bend 12500 -setnumber lc_fault 500 -o sharp015vs_h500.msh
```

| Mesh | Bend | Nodes | Triangles | sha256 |
|---|---|---|---|---|
| straight_h500 | none | 15,620 | 30,854 | cc57a4cae781e1bed98da25f3b4e951624f5135e4b0c6028d6a4573cb63c7f04 |
| sharp015_h500 | 1.5 deg at x = 0 | 15,628 | 30,870 | fb5b5526d68d74106b5146be086e55089aad18e03edcc1d5a7ded6be3d39319d |
| sharp05_h500 | 5 deg at x = 0 | 15,609 | 30,832 | 4f3ab900ad8a46be406c270c69610d588a8c4780e7a35d2dc48b7310e4b269fd |
| arc05L10_h500 | 5 deg over s = -5 to 5 km | 15,644 | 30,902 | 05956834b5a5d81e72192ac56fb83ec1b02079f0ffbbe45eab33d7c92337a529 |
| sharp015vs_h500 | 1.5 deg at x = 12.5 km | 15,647 | 30,908 | 83c995f896981a81060594d26e7081841aa8e8b5158aaca6fbd8d3cfe7ae9d04 |

Each has 60 fault edges of 500 m.

## Length and cost (estimates, not measured)

- **Length:** 12,000 yr, so about 3 events per run. This is an estimate:
  - an earlier Mac run on the straight fault in the +-240 km box had its
    first event at 827 yr and none more by 1,500 yr;
  - the box loads the fault through its sides, so the loading rate falls
    about as 1 / box size; an earlier +-60 km run had its second event at
    1,186 yr, which suggests about 4,000 yr between events here.
- **Wall time:** about 15 to 20 min per case on 4 ranks, about 5 min for a
  run that stops at tension, about 2.5 to 3 h for all twelve (scaled from
  that earlier run: 5.8 min for 1,500 yr with one event).
- **Cap:** `WALL_LIMIT` (default 2700 s) per run. Past it the driver stops
  with a checkpoint (exit 8), and the run can go on later with `--restart`.
- **Disk:** about 0.3 GB per fully dynamic run, much less for QD only.

## Files

- `mesh/*.msh`: the five meshes (above).
- `configs/<case>.toml`: the twelve decks.
- `scripts/make_decks.py`: makes the decks; `--check` compares.
- `scripts/run_simple_mac.sh`: the runs on the Mac, one after the other.
  - Checks: the GitHub-first guard; the binary's code equals HEAD
    (`BUILD_COMMIT_MAC_FDCYCLE`); the decks equal make_decks.py's; the mesh
    sha256s.
  - Takes the Mac run lock; OMP_NUM_THREADS=1; 4 ranks by default.
  - A run that stops at sigma_n' < 0 is a result; the script goes on.
  - Keeps the first 200 `[FRIC-GUARD] sigma_n_eff <= 0` lines in run.log
    and counts the rest.
- `scripts/analyse_simple.py`: `summary.txt`, `summary.json`,
  `simple_cases.pdf`. Per run:
  - end state: exit code, wall time, abort or STOP line, final time;
  - events: time, duration, rupture length, mean slip, peak V;
  - each hand-back: max V dynamic, QD balance, the fault-wide offset
    median(tau_qd - tau_dyn), the largest |tau_qd - tau_dyn| and where;
  - at the bend stations: sigma_n' at t = 0, its lowest value (in an event
    or between), the first time it is <= 0, its value after each event,
    and the slip split into event slip and creep;
  - where the friction law first saw sigma_n' <= 0 (FRIC-GUARD);
  - the lowest sigma_n' on the whole fault (fault snapshots, needs h5py).

## Workflow (Mac, from `miniapps/seas` of the checkout at the pushed commit)

```bash
conda activate mumps-petsc
make -f tests/2d_inplane/2d_inplane.mk seas_fd_cycle_driver
git rev-parse HEAD > BUILD_COMMIT_MAC_FDCYCLE
R=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/simple_cases
simple_cases/scripts/run_simple_mac.sh $R 4
python3 simple_cases/scripts/analyse_simple.py $R
```

## How to read the results (from the plan)

| Compare | If this happens | It means |
|---|---|---|
| Q1, Q2 against Q0 | tension at the bend in QD only | bends cause it from slip alone |
| Q1, Q2 against Q0 | no tension in QD | the dynamic phase is needed, or the slip was too small |
| D0 | hand-back mismatch on the straight fault | a switching problem, not a bend problem |
| D1f, D2f against D0 | mismatch only with bends | the hand-back problem comes from the bend |
| D0, D3, D1f, D2f | all hand back cleanly | the UCERF failure needs something these runs lack (a long fault, large slip) |
| Q2, D2 against Q3, D3 | the arc passes, the sharp bend fails | the corner matters, not the turn |
| D2 against Q2 | similar MPa per metre of slip | the tension follows the slip, whatever the regime |
| Q4, D4 | tension at the creeping bend | creep alone can take a bend into tension |

## Checks done before the runs (cloud, 2026-10-31)

- **Meshes:** made from the repo's `bend.geo`, byte for byte the drafts
  shown in the plan. The nine bend convergence meshes, made again from the
  updated `bend.geo`, are byte for byte the committed ones.
- **Stations:** every station lies on a fault segment of its mesh (largest
  distance off the segment 3e-13 m); `bend_w` and `bend_e` each lie on one
  face only.
- **Start-up:** each of the twelve decks ran 5 QD steps on 2 ranks in the
  cloud (output thrown away; not a result). The binary was the cloud build
  of a4a56fcf, whose fd_cycle code equals HEAD's.
  - All exit 0.
  - Q decks start the QD regime only (P1), D decks the switched one (P4).
  - D1f and D2f print the 1 MPa floor.
  - Every station is found on a fault face at distance 0, and `bend_w` and
    `bend_e` on different faces.
- **`run_simple_mac.sh`** (the loop, with a stub mpirun and without the git
  guard): it records the exit code through the log filter, keeps 200
  FRIC-GUARD lines and counts the rest, writes RUN_INFO.txt, skips an
  existing folder, refuses a held lock and releases its own.
- **`analyse_simple.py`:** on a synthetic run with known answers it
  returned the event slip (4.000 m) and creep (0.110 m), the first
  sigma_n' <= 0 (in event 2), the FRIC-GUARD count (205), and the
  hand-back offset (-0.005 MPa). On the twelve start-up folders it reads
  every station.
