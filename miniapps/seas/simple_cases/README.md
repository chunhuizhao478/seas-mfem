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
questions, so they go first), then D5f, added after D1f, and D5f on two
finer fault meshes (last rows).

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
| D5f_arc015L10_floor | arc015L10_h500 | fully dynamic | 1 MPa | D1f with the bend rounded (added later, below) |
| D5f_arc015L10_h250_floor | arc015L10_h250 | fully dynamic | 1 MPa | D5f with 250 m on the fault (below) |
| D5f_arc015L10_h125_floor | arc015L10_h125 | fully dynamic | 1 MPa | D5f with 125 m on the fault (below) |

- **No floor (all but D1f, D2f):** a run that reaches sigma_n' <= 0 stops
  there, and that time and place is a result.
- **Bent runs without a floor (D1, D2):** they will likely stop in event 1,
  before any hand-back (estimate). D1f and D2f test the hand-back on a
  bend instead.

### D5f: D1f with the bend rounded

Chunhui Zhao, 2026-10-31: "let's prepare for this case, a smooth design
around this bend, and rerun". D1f found sigma_n' at its 1.5 deg corner
growing about 100 MPa per event (-71, -168, -267 MPa on the west face;
+124, +225, +328 MPa on the east face), all at the vertex.

- **Only the corner changes.** The same 1.5 deg turn is spread over a 10 km
  arc centred at x = 0 (x from -5 to +5 km, radius 382 km), as the 5 deg
  arc of Q3/D3. Fault, friction, the 1 MPa floor, box, walls, output and
  stations are D1f's.
- **On the 500 m mesh** the arc is 20 straight segments between points on
  the curve: a turn of 0.075 deg at each of 19 vertices and 0.0375 deg at
  the two arc ends, instead of 1.5 deg at one vertex.
- **What it tests:** if the spike comes from the corner, it should fall to
  a few MPa per event, spread along the arc, with no tension and no floor
  used (estimate).

### D5f refined: 250 m and 125 m on the fault

Chunhui Zhao, 2026-10-31: "okay keep it recorded, and let's refine the mesh
1 level / 2 levels for fault and rerun this problem". D5f's arc is a polygon
of 500 m straight elements, and sigma_n' zigzags between about +7 and +45
MPa from vertex to vertex (Results, below).

- **Only lc_fault changes** (500 to 250 and 125 m). The points stay on the
  true arc, so the turn at a vertex halves each level: 0.0375 deg at 250 m,
  0.019 deg at 125 m. The mesh grows from the fault as in `bend.geo`; the
  box, the far field (5 km) and everything in the deck stay D5f's.
- **What it answers:** does the vertex-to-vertex swing shrink with the
  turn (about half per level if each vertex acts as a small corner), and
  does the west tip's tension converge or grow as h shrinks?
- **Stations:** the same places; `bend_w` and `bend_e` follow each mesh's
  chord, 1 m from the middle vertex.

## Setup

All the decks are `config/fd_cycle/mapview_rs_fd.toml` with the edits
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
- **Floor (D1f, D2f, D5f):** `[cycle] sigma_n_floor_pa = 1.0e6` (UCERF run8's).

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
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 1.5 -setnumber arc 10000 -setnumber s_bend 0     -setnumber lc_fault 500 -o arc015L10_h500.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 1.5 -setnumber arc 10000 -setnumber s_bend 0     -setnumber lc_fault 250 -o arc015L10_h250.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 1.5 -setnumber arc 10000 -setnumber s_bend 0     -setnumber lc_fault 125 -o arc015L10_h125.msh
```

| Mesh | Bend | Nodes | Triangles | sha256 |
|---|---|---|---|---|
| straight_h500 | none | 15,620 | 30,854 | cc57a4cae781e1bed98da25f3b4e951624f5135e4b0c6028d6a4573cb63c7f04 |
| sharp015_h500 | 1.5 deg at x = 0 | 15,628 | 30,870 | fb5b5526d68d74106b5146be086e55089aad18e03edcc1d5a7ded6be3d39319d |
| sharp05_h500 | 5 deg at x = 0 | 15,609 | 30,832 | 4f3ab900ad8a46be406c270c69610d588a8c4780e7a35d2dc48b7310e4b269fd |
| arc05L10_h500 | 5 deg over s = -5 to 5 km | 15,644 | 30,902 | 05956834b5a5d81e72192ac56fb83ec1b02079f0ffbbe45eab33d7c92337a529 |
| sharp015vs_h500 | 1.5 deg at x = 12.5 km | 15,647 | 30,908 | 83c995f896981a81060594d26e7081841aa8e8b5158aaca6fbd8d3cfe7ae9d04 |
| arc015L10_h500 | 1.5 deg over s = -5 to 5 km | 15,637 | 30,888 | ff417d7642cfa57c981dcab948e55f69d4a6ed1883cb64743b11e88d88ecf9ce |
| arc015L10_h250 | the same, 250 m on the fault | 19,982 | 39,578 | 6382f4f3d48496b9783079b4f201f0d66392acf5e0a30c36ebad28fc767d6fce |
| arc015L10_h125 | the same, 125 m on the fault | 27,355 | 54,324 | ad2bcffb83e1a409e63d1173076bdfdbcacb5a7a86777bb6aa39e857183ffe9f |

The 500 m meshes have 60 fault edges; arc015L10_h250 120 and arc015L10_h125
240.

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
- **Measured (Mac, 932ccc20):** D1f took 2,419 s for 3 events about 4,600
  yr apart, D2f 2,339 s; so about 40 min per fully dynamic case, not 15 to
  20. D5f is run with `WALL_LIMIT=3600` so the cap does not cut it short.
- **D5f's split (measured, run.log):** 2,498 s in all; the 3 dynamic
  phases about 96 s (about 23,000 steps each at dt 3.8e-3 s, 1.4 ms per
  step); QD about 2,400 s (775,592 RHS stages, about 3.1 ms each, with
  output). G (N_f = 120) builds in 0.3 s.
- **Refined D5f (estimate, not measured):** the dynamic dt halves each
  level and the mesh grows 1.28x and 1.75x, so the dynamic phases take about
  4 and 11 min. QD has 2x and 4x the fault nodes (N_f 240 and 480); if its
  stage cost and its stage count each grow by up to 1.5x per level, QD
  takes about 40 to 90 min at 250 m and 40 to 200 min at 125 m. So about 1
  to 1.5 h and 1 to 3.5 h. They run with `WALL_LIMIT=21600` (6 h each), one
  after the other.

## Files

- `mesh/*.msh`: the eight meshes (above).
- `configs/<case>.toml`: the fifteen decks.
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
make -f tests/2d_inplane/2d_inplane.mk SEAS_FS_LIB= seas_fd_cycle_driver   # macOS: no -lstdc++fs
git rev-parse HEAD > BUILD_COMMIT_MAC_FDCYCLE
R=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/simple_cases
simple_cases/scripts/run_simple_mac.sh $R 4
WALL_LIMIT=3600 simple_cases/scripts/run_simple_mac.sh $R 4 D5f_arc015L10_floor   # D5f alone
WALL_LIMIT=21600 simple_cases/scripts/run_simple_mac.sh $R 4 D5f_arc015L10_h250_floor D5f_arc015L10_h125_floor
python3 simple_cases/scripts/analyse_simple.py $R
```

## Results (Mac, 2026-10-31)

Run with `run_simple_mac.sh` on 4 ranks; folders
`/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/simple_cases/<case>_<commit>/`;
summary `compare/summary.txt`. sigma_n' in MPa (compression positive,
prestress 25); hand-back = QD balance max V (m/s) / the fault-wide offset
median(tau_qd - tau_dyn) (MPa).

| Case | Commit | Exit, wall | Events (yr) | Hand-backs |
|---|---|---|---|---|
| D1f | 932ccc20 | 0, 2,419 s | 843, 5,391, 10,063 (7.3, 7.0, 7.2 m) | 8.4e-4 / -0.043; 8.4e-4 / -0.042; 8.5e-4 / -0.044 |
| D2f | 932ccc20 | 1, 2,339 s | 854, 6,007 (6.7, 7.4 m) | 8.5e-4 / -0.060; 8.6e-4 / -0.051 |
| Q1 | 932ccc20 | 8 (cap), 2,702 s | event 1 only | QD only |
| D5f | fc17fbd3 | 0, 2,498 s | 846, 5,345, 9,930 (7.5, 7.5, 7.6 m) | 8.3e-4 / -0.044; 8.3e-4 / -0.045; 8.4e-4 / -0.045 |

Q0, Q2, Q3, Q4 and D0 to D4 were not run (cancelled 2026-10-31: "cancel all
other queued job").

- **Hand-backs are clean in every run,** bent or not (UCERF's: 0.253 m/s
  and -0.78 MPa). A bend alone does not reproduce the UCERF hand-back.
- **D1f, the sharp 1.5 deg corner:** sigma_n' min / max per event, from
  every fault snapshot, are -71 / +124, -168 / +225, -267 / +328 MPa,
  both on the two nodes at the vertex (min on the west face, max on the
  east), reached 15 to 22 s into each phase. They grow about 100 MPa per
  event and fall back to about 25 MPa within about 1 km.
- **D2f, the sharp 5 deg corner:** the vertex was at -517 / +608 MPa after
  2 events. At about 11,749 yr, in calm QD (dt about 5.5e4 s), it stopped
  with `Non-finite V from friction solver at DOF 67`
  (rate_state_fault.hpp:1378) at a floored node 0.5 km from the corner:
  tau 6.4e10 Pa and psi -32 in a trial stage. Friction there is zero, the
  root is at the bracket's upper end and rounding fails both brackets
  (dieterich_ruina.hpp). The runaway stage before it is the open question;
  its checkpoint (step 130,000) is 475 steps before the abort.
- **Q1, QD only:** the corner goes tensile in event 1 (bend_w -50.6 MPa),
  so the tension needs no dynamic phase. With no floor the tensile node
  kept slipping at about 1e-5 m/s, and QD stalled at dt about 0.04 s at
  843.2 yr until the cap.
- **D5f, the 10 km arc:** no tension at the bend; its min / max per event
  are +18 / +32, +12 / +39, +6 / +47 MPa (7x below D1f). But sigma_n'
  zigzags between about +7 and +45 MPa from vertex to vertex along the
  arc, and the swing still grows each event. The west tip (DOF 22, x = -15
  km) goes tensile at the end of each phase: +17, -0.9, -18.2 MPa, the
  only node that used the floor (D1f's west tip: +18 to +39 MPa; the east
  tip stays near 25 MPa). Not yet explained.
- **Figures (Mac):** `<case>_<commit>/gifs/` (fault and domain GIFs,
  static min / max figures); pages "D1f GIFs" and "D1f vs D5f".

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
| D5f against D1f | the arc's sigma_n' change is a few MPa per event, no floor used | the spike is the corner's, not the 1.5 deg turn's |

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

### D5f (cloud, 2026-10-31)

- **Mesh:** the commands above give the committed `sharp015_h500.msh` byte
  for byte again; `arc015L10_h500.msh` comes from the same `bend.geo`. Its
  61 fault points lie on the arc within 4e-12 m; the turns sum to 1.5 deg,
  at most 0.075 deg at a vertex, all within x = -5 to +5 km.
- **Stations:** every fault station of D5f lies on a fault segment (largest
  distance 6e-14 m).
- **Decks:** `make_decks.py --check` gives "same" for all thirteen; the
  twelve earlier decks are unchanged.
- **No start-up run in the cloud** (Chunhui Zhao, 2026-10-30: all runs on
  the Mac). D5f's first steps on the Mac are its start-up check.

### D5f refined (cloud, 2026-10-31)

- **Meshes:** from the same `bend.geo` and gmsh 4.15.2; 121 and 241 fault
  points on the arc; node and triangle counts above.
- **Stations:** every fault station of the two decks lies on a fault
  segment (at most 6e-14 m off).
- **Decks:** `make_decks.py --check` gives "same" for all fifteen; the
  thirteen earlier decks are unchanged. The two new decks differ from
  D5f's only in the header, the mesh path, the prefix, two comments
  (fault size, N_f) and `bend_w` / `bend_e` (each mesh's chord).
- **No start-up run in the cloud.**
