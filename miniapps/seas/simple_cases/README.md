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
questions, so they go first), then D5f, added after D1f, D5f on two
finer fault meshes, the three P0 penalty runs, the five runs of the
vertex-normal plan's P3 and P4, the five smoothed bends, and the three
tip zones (last rows).

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
| D5f_arc015L10_pf1_ev1 | arc015L10_h500 | fully dynamic, event 1 | 1 MPa | P0: penalty_factor 1 at cfl 0.25 (below) |
| D5f_arc015L10_pf2_ev1 | arc015L10_h500 | fully dynamic, event 1 | 1 MPa | P0: penalty_factor 2 (below) |
| D5f_arc015L10_pf4_ev1 | arc015L10_h500 | fully dynamic, event 1 | 1 MPa | P0: penalty_factor 4 (below) |
| D0_straight_ev1 | straight_h500 | fully dynamic, event 1 | none | P3: the straight control, face normals (below) |
| D0_straight_avg_ev1 | straight_h500 | fully dynamic, event 1 | none | P3: the same, averaged vertex normals |
| D5f_arc015L10_avg_ev1 | arc015L10_h500 | fully dynamic, event 1 | 1 MPa | P3: D5f's event 1, averaged (smoke) |
| D5f_arc015L10_avg_floor | arc015L10_h500 | fully dynamic | 1 MPa | P4: D5f, averaged |
| D5f_arc015L10_h250_avg_floor | arc015L10_h250 | fully dynamic | 1 MPa | P4: D5f at 250 m, averaged |
| D7f_arc08L10_h250_avg_floor | arc08L10_h250 | fully dynamic | 1 MPa | smoothed bend: 8 deg over 10 km, averaged (below) |
| D9f_arc08L2p5_h250_avg_floor | arc08L2p5_h250 | fully dynamic | 1 MPa | smoothed bend: 8 deg over 2.5 km, averaged |
| D6f_arc03L10_h250_avg_floor | arc03L10_h250 | fully dynamic | 1 MPa | smoothed bend: 3 deg over 10 km, averaged |
| D8f_arc03L2p5_h250_avg_floor | arc03L2p5_h250 | fully dynamic | 1 MPa | smoothed bend: 3 deg over 2.5 km, averaged |
| D10f_sbend03L05_h250_avg_floor | sbend03L05_h250 | fully dynamic | 1 MPa | smoothed bends: +3 then -3 deg, 5 km arcs, averaged |
| D5f_arc015L10_h250_avg_floor_T1a08 | arc015L10_h250 | fully dynamic | 1 MPa | tip zone T1: a = 0.08 within 3 km of each tip (below) |
| D5f_arc015L10_h250_avg_floor_T2w8 | arc015L10_h250 | fully dynamic | 1 MPa | tip zone T2: a = 0.04 within 6 km of each tip |
| D5f_arc015L10_h250_avg_floor_T3a08w8 | arc015L10_h250 | fully dynamic | 1 MPa | tip zone T3: a = 0.08 within 6 km of each tip |

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

### P0: the penalty scan (vertex-normal plan)

Chunhui Zhao, 2026-11-01: "send this document to 2d thread and start
implementation and test", for
`/mnt/project-files/2d-inplane/nonplanar-literature/vertex-normal-plan-rev1.pdf`
(averaged vertex normals on the fault). Its phase P0 checks the mechanism
before any code changes.

- **The inferred mechanism:** at an arc vertex the two straight faces
  prescribe the slip along two tangents a turn apart, and the Nitsche
  penalty acts on the mismatch. If so, the vertex jump in sigma_n'
  (about 13 MPa in event 1 at 500 m) grows with `penalty_factor`.
- **Only three keys change** from D5f's deck:
  - `[cycle] penalty_factor` = 1, 2, 4;
  - `[numerics] cfl` = 0.25 (D5f: 0.5). The explicit step does not see the
    penalty (`solver/nitsche_cfl.hpp`), but the face stiffness does, so
    omega_max can grow up to about 2x at factor 4 (estimate). 0.25 keeps
    D5f's margin at factor 4;
  - `[cycle.switch] max_dynamic_phases = 1`: the run stops at the end of
    event 1, before switching back to QD.
- **The factor 1 run** uses the same cfl 0.25, so the step's own effect is
  separated from the penalty's (D5f itself is factor 1 at cfl 0.5).
- **What it answers** (analyse_simple.py, "Vertex jumps"):
  - the jump grows with the factor (about in proportion, if the penalty
    sets it): the mechanism holds, and P1 starts;
  - the jump does not change: stop and look at the mechanism again.

### P3 and P4: averaged vertex normals (vertex-normal plan)

The same request as P0. P1 to P3 add `[fault_geometry] vertex_normals =
"averaged"` (default `"face"`): at a fault vertex where two faces turn by
at most 5 deg, both faces use n_v = (n_A + n_B) / |n_A + n_B|, and the
frame varies along each face between its two end normals. Each `_avg`
deck is its base deck plus that one line.

- **P3, regression and smoke** (each to the end of event 1,
  `max_dynamic_phases = 1`):
  - `D0_straight_ev1` against D0's run (38e7ad28): the option off is
    unchanged (the plan asks for bit for bit; on more than one rank two
    runs of one build can differ in the last bits, so equal to round-off);
  - `D0_straight_avg_ev1` against `D0_straight_ev1`: on a straight fault
    the option changes nothing (no face frame varies, so both modes take
    the same code path; on 4 ranks two runs of one deck differ by about
    1e-7 after 827 yr, so the bit-for-bit check is on 1 rank, below);
  - `D5f_arc015L10_avg_ev1`: the run is clean and the event-1 vertex jump
    is far below D5f's 13.1 MPa (the target: under 1 MPa).
- **P4** (12,000 yr, about 3 events, as D5f): `D5f_arc015L10_avg_floor`
  and `D5f_arc015L10_h250_avg_floor` against D5f and D5f 250 m. The
  plan's targets: the vertex jump under 1 MPa per event; event times
  within 1% and mean slip within 2% of face mode; hand-backs as clean; no
  floored DOF.
- **Output:** with the option on the driver also writes
  `<prefix>_vertices.txt`, one row per fault vertex: key, x, y, kind (0
  tip, 1 smooth, 2 kink), turn angle (deg) and n_v.

### Smoothed bends: how tight can a smoothed bend be?

Chunhui Zhao, 2026-11-01, after P4: "do you suggest we run a few more
cases by adjusting the angles (similar to what UCERF fault) before we move
back to UCERF fault configuration?"; "we can smooth the trace bends if
needed, that's not an issue I think".

- **Question:** the UCERF traces will be smoothed, so what arc is tight
  enough to stay clean with averaged normals? The UCERF trace at 500 m
  has about 180 vertices turning 1 to 5 deg and 5 above 5 deg (the
  largest 8.74 deg).
- **Cases:** P4's 250 m deck (averaged, 1 MPa floor, 12,000 yr) on five
  meshes; each arc's ends and middle are fault points.

| Case | Bend | Radius | Turn per vertex |
|---|---|---|---|
| D6f_arc03L10_h250_avg_floor | 3 deg over s = -5 to 5 km | 191 km | 0.075 deg |
| D7f_arc08L10_h250_avg_floor | 8 deg over s = -5 to 5 km | 71.6 km | 0.2 deg |
| D8f_arc03L2p5_h250_avg_floor | 3 deg over s = -1.25 to 1.25 km | 47.7 km | 0.3 deg |
| D9f_arc08L2p5_h250_avg_floor | 8 deg over s = -1.25 to 1.25 km | 17.9 km | 0.8 deg |
| D10f_sbend03L05_h250_avg_floor | +3 deg over s = -6.25 to -1.25 km, -3 deg over 1.25 to 6.25 km | 95.5 km | 0.15 deg |

- **Comparisons:** D6f against D8f and D7f against D9f change the arc
  length at one angle; D6f against D7f the angle at one length; D10f two
  bends in a row (the fault ends parallel to its start, 0.39 km north).
  All against P4's 1.5 deg arc at 250 m (`D5f_arc015L10_h250_avg_floor`).
- **What to read:** the vertex jump (|x| <= 5 km; 7 km for D10f), the arc's
  sigma_n' min and max and how they change per event, any tension or
  floored DOF, event times and mean slip, the hand-backs.
- **Expected (estimate, not tested):** the smooth drift across the arc
  (about 1.5 MPa per event at each end of P4's 1.5 deg arc) is the bend's
  own stress and grows with the angle; it scales with slip times
  curvature, so the 17.9 km arc may reach tension in a few events (about
  13 MPa per event, rough) with clean numerics. That would set the
  smoothing limit; it is physics, not mesh.
- **Not covered:** UCERF events slip much more than the 7.5 m here, and
  the bend stress grows with slip; the UCERF hand-back problem does not
  appear in any simple case so far.
- **Order:** the 8 deg cases first (D7f, D9f), then D6f, D8f, D10f.

### Tip zones: a velocity-strengthening patch measured from each tip

Chunhui Zhao, 2026-11-01: "I would like to add functionality about
velocity strengthening patch near the tip of the fault, and test on
previous 1.5 degree fault with averaged normal, and see if this could
reduce the spikes of normal stress on both fault ends".

- **The functionality:** `[cycle.fields.a] kind = "tip_zone"`
  (`cycle/fault_chains.hpp`, 2-D only), with d the distance along the
  DOF's own fault to its nearer free tip:

  | d | a |
  |---|---|
  | d <= width - taper | tip |
  | width - taper < d < width | linear from tip to inside |
  | d >= width | inside |

  - `max_fraction` (optional) caps the width at max_fraction x L per fault.
  - A fault end on the outer boundary is not a tip.
  - Junctions and closed faults are refused.
  - The driver logs the faults and tips it found.
  - For UCERF's many curved faults, where a rectangle in x cannot follow
    each fault.
- **Already in every simple-case deck:** the rect taper gives a = 0.04
  (a - b = +0.01) for |x| >= 12 km, from 0.004 at 10 km. On this 30 km
  fault that is the tip zone width 5 km, taper 2 km, tip 0.04. P4's 250 m
  run still drifts at both ends: by event 3, 12.1 MPa at x = 14.5 km and
  43.8 MPa at x = -14.75 km (from 25), 7 to 8 MPa per event. The straight
  D0 drifts the same way (12.3 / 30.4 MPa).
- **Cases:** P4's 250 m deck (`D5f_arc015L10_h250_avg_floor`: averaged,
  1 MPa floor, 12,000 yr) with the tip zone in place of the rect taper:

  | Case | tip a (a - b) | Full tip value | Ramp | Weakening part |
  |---|---|---|---|---|
  | T1 `_T1a08` | 0.08 (+0.05) | d <= 3 km | 3 to 5 km | 20 km |
  | T2 `_T2w8` | 0.04 (+0.01) | d <= 6 km | 6 to 8 km | 14 km |
  | T3 `_T3a08w8` | 0.08 (+0.05) | d <= 6 km | 6 to 8 km | 14 km |

  - b = 0.03 everywhere.
  - The prestress uses the local a ([cycle.prestress] a = "local"), so
    the zone starts at steady state at V_p for its own a (tau / sigma_n'
    0.255 at a = 0.08, 0.531 at 0.04).
- **What to read:** the tip ends per dynamic phase (analyse_simple.py):
  sigma_n' at each tip DOF and the largest change within 1 km of it, the
  tip's slip and peak V. Also event times, mean slip, the hand-backs and
  the arc's vertex jump. All against `D5f_arc015L10_h250_avg_floor`.
- **Expected (estimate, not tested):**
  - The drift should drop with the slip near the tips.
  - It will not vanish, since the fault still ends in intact rock and
    creeps there between events.
  - T2 and T3 shrink the weakening part from 20 to 14 km (about 2
    nucleation lengths), so their events change too.
- **Order:** T1, T2, T3, after the smoothed bends.

### Relaxation: D9f with sigma_n' relaxation and a 10 MPa floor

Chunhui Zhao, 2026-11-02: "implement the relaxation, test and verify it,
and apply to UCERF 2D plane strain, with 10MPa cap, and rerun the problem".
Design: `/mnt/project-files/2d-inplane/relaxation/relaxation-design-rev0.md`
(Ozawa, Ando & Dunham 2023, eq. 6, normal stress only).

- **The functionality:** `[cycle.relaxation] t_relax_yr`
  (`cycle/normal_stress_relaxation.hpp`), off when absent:

  ```
  sigma_law = max(sigma0 + dsig - r, sigma_floor)      both regimes
  dr/dt     = (dsig - r) / t_r                        QD regime; r = 0 at t = 0
  ```

  - dsig: the elastic normal-traction change (slip plus loading).
  - r is frozen in the dynamic regime (the normal-stress hook).
  - The bulk and the outputs keep the elastic stress; the fault
    snapshots add `normal_stress_relaxation` (r, MPa).
- **Cases:** D9f (8 deg over 2.5 km, R 17.9 km, the tightest smoothed
  bend; 2 events in 12,000 yr, 3 DOFs on the 1 MPa floor at event 2),
  with the floor at 10 MPa (the UCERF rerun's) and 30,000 yr (about 5
  events, estimate):

  | Case | Floor | Relaxation |
  |---|---|---|
  | `D9f_arc08L2p5_h250_avg_floor10` | 10 MPa | none |
  | `D9f_arc08L2p5_h250_avg_floor10_relax5k` | 10 MPa | t_r = 5,000 yr |

- **What to read:** sigma_n' of the law (`normal_stress` -
  `normal_stress_relaxation`) at the arc ends per event: with relaxation
  it should level off, without it keep falling to the floor; on DOFs
  away from the bend, r against a (t - t_r (1 - exp(-t / t_r))) of the
  load drift a. Also event times and slip, the floored DOFs, the
  hand-backs.
- **Expected (estimate, not tested):** the D9f arc-end notch (about
  15-30 MPa per event) is large next to 25 MPa, so the 10 MPa floor still
  acts in both runs; relaxation should shorten how long DOFs sit on it.

## PML hand-back check (PML design plan rev2, P3, test T5; 2026-11-04)

- **Why:** the plan's Proposition 3: with the CFS-PML the dynamic phase
  settles to the QD static state, so the hand-back jump
  |tau_qd - tau_dyn| should fall to the D4 target, 0.05 MPa.
- **Cases:** D0 (straight), event 1 (about 827 yr) and its hand-back,
  900 yr; Kelvin-Voigt 0.1 as D0. The band: every outer face, 150 km
  (30 elements of 5 km), N 2, R0 1e-3.

  | Case | Outer boundary in the dynamic regime |
  |---|---|
  | `D0_straight_hb_ref` | absorbing walls (D0's) |
  | `D0_straight_hb_pml0p005` | CFS-PML, f_alpha 0.005 Hz (s0 14.2, dt / 3.8) |
  | `D0_straight_hb_pml0p01` | CFS-PML, f_alpha 0.01 Hz (s0 7.6, dt / 2.8) |
  | `D0_straight_hb_pml0p02` | CFS-PML, f_alpha 0.02 Hz (s0 4.3, dt / 2.1) |
  | `D0_straight_hb_sponge` | sponge comparator (same profile, no stretch) |

- **What to read:** the hand-back line ("largest |tau_qd - tau_dyn|") and
  `<prefix>_handback_phase1.txt`; the QD balance V; event time and slip
  against the reference; the dynamic phase's length and wall time.
- **Expected (estimate):** the CFS-PML cases below the reference's jump,
  the smallest f_alpha lowest; the sponge between.

## PML on D9f (PML design plan rev2, P4, test T6; 2026-11-04)

- **Why:** D9f's hand-backs carry the wall offset (about -0.05 MPa,
  -(1.1-1.7) mu P / A); with Dirichlet walls and the band in both regimes
  it should go (Proposition 3).
- **Case:** `D9f_arc08L2p5_h250_avg_floor10_relax5k_pml0p01` = relaxation
  V1 (t_r 5 kyr, 10 MPa floor) with the CFS-PML band (150 km, f_alpha
  0.01 Hz, D3 auto), 12,000 yr (its first two events). Reference: the
  relax5k run's first two events (absorbing walls).
- **What to read:** each hand-back's largest |tau_qd - tau_dyn| and the
  fault median of tau_qd - tau_dyn (the wall offset), QD balance V, event
  times and slip, floored DOFs, sigma_n' at the arc ends.

## Fault opening A on D9f (fault opening design rev 1.1, test T5; 2026-11-04)

- **Why:** the 10 MPa floor control (`D9f_arc08L2p5_h250_avg_floor10`)
  ratchets the elastic sigma_n' at the arc ends to tension event by event;
  the floor then holds strength there that a tensile point does not have.
- **Case:** `D9f_arc08L2p5_h250_avg_floor10_open` = that control with
  `[cycle.opening] mode = "shear_free"` (version A, cycle/fault_opening.hpp):
  a DOF whose elastic sigma_n' is tensile is open, shear-free in both
  regimes, psi frozen, left out of max V and the hand-back guard; the faces
  stay glued; 30,000 yr. Reference: the control's run.
- **What to read:** the open DOFs per hand-back (driver log, hand-back
  table columns open_dyn and open_qd), exit code, events and their slip,
  the hand-back jump at closed DOFs, the elastic sigma_n' at the arc ends
  (does it still ratchet), and whether the open set chatters (settle
  flips per accepted step).

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
(its new inputs `box` and `s_bend`, and for the smoothed bends' pair
`theta2_deg`, `arc2` and `s_bend2`; its defaults still give the nine bend
meshes byte for byte, and all thirteen meshes here were remade byte for
byte after the second bend was added):

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
# the smoothed bends (2026-11-01)
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 3   -setnumber arc 10000 -setnumber s_bend 0     -setnumber lc_fault 250 -o arc03L10_h250.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 8   -setnumber arc 10000 -setnumber s_bend 0     -setnumber lc_fault 250 -o arc08L10_h250.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 3   -setnumber arc 2500  -setnumber s_bend 0     -setnumber lc_fault 250 -o arc03L2p5_h250.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 8   -setnumber arc 2500  -setnumber s_bend 0     -setnumber lc_fault 250 -o arc08L2p5_h250.msh
gmsh -2 $G -format msh22 -setnumber box 240e3 -setnumber theta_deg 3   -setnumber arc 5000  -setnumber s_bend -3750 -setnumber theta2_deg -3 -setnumber arc2 5000 -setnumber s_bend2 3750 -setnumber lc_fault 250 -o sbend03L05_h250.msh
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
| arc03L10_h250 | 3 deg over s = -5 to 5 km, 250 m | 19,978 | 39,570 | f535f8fecae6fdcf5265d18e5d442ff3ec13b1d2949a9ca47f6cdbb44882bdd1 |
| arc08L10_h250 | 8 deg over s = -5 to 5 km, 250 m | 20,004 | 39,622 | 102383a4cbfd420a3c21cfdddafe9e6b5ef37e49e9bbeb121b55c32336eb70e9 |
| arc03L2p5_h250 | 3 deg over s = -1.25 to 1.25 km, 250 m | 19,957 | 39,528 | f497007515d145c22cb67b67fa0f4756fff5777bca00b0a22e5eca80d106bda5 |
| arc08L2p5_h250 | 8 deg over s = -1.25 to 1.25 km, 250 m | 19,999 | 39,612 | 9a182625e106083f2331b96a26600fe8f39b0d6da8edd907df99854527078505 |
| sbend03L05_h250 | +3 then -3 deg, 5 km arcs, 250 m | 19,979 | 39,572 | 85e6ff91e7d6d8d73664c6f0f198aa965238f0f40eb509d0ceb1eca1d147b376 |

The 500 m meshes have 60 fault edges; the 250 m meshes 120 and
arc015L10_h125 240. On the smoothed-bend meshes every fault point lies on
its curve to 6e-12 m, the turns add up to 3, 8 and 0 deg, and every
station lies on a fault edge to 3e-13 m (checked in the cloud).

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
- **P0 runs (estimate, not measured):** QD to event 1 at about 846 yr is
  about a third of D5f's QD stages (about 13 min); the dynamic phase at
  cfl 0.25 is about 47,000 steps (about 1 min). So about 15 min per run,
  45 min for the three, with `WALL_LIMIT=3600`.
- **P3 and P4 runs (estimate, not measured):** the `_ev1` runs as P0's
  at cfl 0.5, about 15 min each; D5f averaged as D5f, about 42 min; at
  250 m as D5f 250 m, about 53 min (3,189 s).
- **Smoothed bends (estimate, not measured):** as P4's 250 m run (3,120 s
  measured), about 50 min each, about 4.2 h for the five, with
  `WALL_LIMIT=21600`.
- **Tip zones (estimate, not measured):** the same, about 50 min each,
  about 2.5 h for the three, with `WALL_LIMIT=21600`.
- **Refined D5f (estimate, not measured):** the dynamic dt halves each
  level and the mesh grows 1.28x and 1.75x, so the dynamic phases take about
  4 and 11 min. QD has 2x and 4x the fault nodes (N_f 240 and 480); if its
  stage cost and its stage count each grow by up to 1.5x per level, QD
  takes about 40 to 90 min at 250 m and 40 to 200 min at 125 m. So about 1
  to 1.5 h and 1 to 3.5 h. They run with `WALL_LIMIT=21600` (6 h each), one
  after the other.

## Files

- `mesh/*.msh`: the thirteen meshes (above).
- `configs/<case>.toml`: the forty decks.
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
  - the lowest sigma_n' on the whole fault (fault snapshots, needs h5py);
  - the vertex jumps per dynamic phase: at each interior fault vertex the
    largest minus the smallest sigma_n' of its DOFs (one per face), the
    largest over the phase for |x| <= 5 km and over the whole fault, and
    the |x| <= 5 km one at the phase's end;
  - the tip ends per dynamic phase: at the west and east tip DOFs,
    sigma_n' minus its t = 0 value at the phase's end, the slip there and
    the peak |V| in the phase; the largest |sigma_n' - sigma_n'(0)| within
    1 km of each tip (where, signed).

## Workflow (Mac, from `miniapps/seas` of the checkout at the pushed commit)

```bash
conda activate mumps-petsc
make -f tests/2d_inplane/2d_inplane.mk SEAS_FS_LIB= seas_fd_cycle_driver   # macOS: no -lstdc++fs
git rev-parse HEAD > BUILD_COMMIT_MAC_FDCYCLE
R=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/simple_cases
simple_cases/scripts/run_simple_mac.sh $R 4
WALL_LIMIT=3600 simple_cases/scripts/run_simple_mac.sh $R 4 D5f_arc015L10_floor   # D5f alone
WALL_LIMIT=21600 simple_cases/scripts/run_simple_mac.sh $R 4 D5f_arc015L10_h250_floor D5f_arc015L10_h125_floor
WALL_LIMIT=3600 simple_cases/scripts/run_simple_mac.sh $R 4 D5f_arc015L10_pf1_ev1 D5f_arc015L10_pf2_ev1 D5f_arc015L10_pf4_ev1   # P0
WALL_LIMIT=3600 simple_cases/scripts/run_simple_mac.sh $R 4 D0_straight_ev1 D0_straight_avg_ev1 D5f_arc015L10_avg_ev1   # P3
WALL_LIMIT=21600 simple_cases/scripts/run_simple_mac.sh $R 4 D5f_arc015L10_avg_floor D5f_arc015L10_h250_avg_floor   # P4
WALL_LIMIT=21600 simple_cases/scripts/run_simple_mac.sh $R 4 D7f_arc08L10_h250_avg_floor D9f_arc08L2p5_h250_avg_floor   # smoothed bends, 8 deg
WALL_LIMIT=21600 simple_cases/scripts/run_simple_mac.sh $R 4 D6f_arc03L10_h250_avg_floor D8f_arc03L2p5_h250_avg_floor D10f_sbend03L05_h250_avg_floor
WALL_LIMIT=21600 simple_cases/scripts/run_simple_mac.sh $R 4 D5f_arc015L10_h250_avg_floor_T1a08 D5f_arc015L10_h250_avg_floor_T2w8 D5f_arc015L10_h250_avg_floor_T3a08w8   # tip zones
python3 simple_cases/scripts/analyse_simple.py $R
```

## Results (Mac, 2026-10-31 and 2026-11-01)

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
| D5f 250 m | 38e7ad28 | 0, 3,189 s | 843, 5,327, 9,898 (7.5, 7.5, 7.7 m) | 8.4e-4; 8.4e-4; 8.4e-4 |
| D5f 125 m | 38e7ad28 | 1 (cancelled), 2,968 s | 841, 5,315 (7.5, 7.5 m) | 8.4e-4; 8.5e-4 |
| D0 | 38e7ad28 | 0, 2,498 s | 827, 5,202, 9,572 (7.5, 7.3, 7.4 m) | 8.3e-4; 8.2e-4; 8.2e-4 |
| Q0 | 38e7ad28 | 0, 2,409 s | 827, 4,258, 7,700, 11,143 (5.7, 5.6, 5.6, 5.6 m) | QD only |

Q2, Q3, Q4 and D1 to D4 were not run (cancelled 2026-10-31: "cancel all
other queued job"). D0 and Q0 ran later (2026-10-31: "can we also run the
straight lines benchmark"); the 125 m run was stopped after event 2
(2026-10-31: "you can then cancel 125m run, because we don't see much
improvement with mesh refinement"). For the four 38e7ad28 rows the hand-back
column is the QD balance only; their offsets are about -0.04 MPa, as above.

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
- **D5f refined, the vertex jump does not shrink with h:** the largest
  sigma_n' jump across a vertex with |x| <= 5 km, over every dynamic
  snapshot, per event: 13.1 / 26.4 / 40.0 MPa at 500 m, 13.6 / 27.6 / 41.9
  at 250 m, 14.1 / 28.4 at 125 m (2 events). Event times, mean slip and
  hand-backs agree within about 1%.
- **The 500 m west-tip tension was the mesh:** at 250 m the west tip's
  lowest sigma_n' is +23.3 / +18.6 / +13.8 MPa, and the floor is never
  used at 250 or 125 m.
- **D0, the straight control:** the same jump is 0.6 / 0.9 / 1.2 MPa, so
  the sawtooth comes from the arc's vertices. Hand-backs are as clean as
  on the bends.
- **The ends drift even on a straight fault:** sigma_n' at x = -14.5 km
  falls from 25 to 12.3 MPa by event 3 in D0 (9.6 MPa by event 4 in Q0),
  with the tip node rising to about 30 MPa; the east end changes less. Not
  yet explained.
- **QD only against fully dynamic (Q0 against D0):** 4 events about 3,440
  yr apart with mean slip 5.6 m, against 3 events about 4,370 yr apart
  with 7.4 m.
- **P3 and P4, averaged vertex normals (83eeb186, 2026-11-01):**
  - P3: `D0_straight_ev1`, `D0_straight_avg_ev1` and
    `D5f_arc015L10_avg_ev1` exit 0 (260, 212, 208 s); event 1 at 827.01,
    827.01 and 845.67 yr. On 1 rank, D0 averaged and face give identical
    data (every text output and snapshot dataset). On 4 ranks they differ
    by about 1e-7 relative, as much as a repeat of the same run.
  - P4, event times and mean slip within 0.07% and 0.2% of face mode,
    hand-backs as clean. The vertex jump per event: 0.44 / 0.73 / 0.99 MPa
    at 500 m (face 13.07 / 26.44 / 40.05), 0.34 / 0.52 / 0.70 at 250 m (face
    13.64 / 27.61 / 41.88). The arc's sigma_n' at 250 m: 24.2-26.0,
    23.4-27.7, 22.6-29.4 MPa per event, the same at 500 m to about 0.3 MPa.
  - The 500 m west tip is floored in events 2 and 3, as in face mode; at
    250 m nothing is floored. Near the ends at 250 m sigma_n' reaches +12
    and +44 MPa by event 3 (as D0's end drift).
- **Figures (Mac):** `<case>_<commit>/gifs/` (fault and domain GIFs,
  static min / max figures); pages "D1f GIFs", "D1f vs D5f", "D5f arc
  refinement" and the summary "Simple Fault Cases"
  (https://claude.ai/artifact/VWN8Gc4QZW8XLYTQvUYdHb). The jumps above are
  the Mac's measure; analyse_simple.py's "Vertex jumps" computes the same.

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
