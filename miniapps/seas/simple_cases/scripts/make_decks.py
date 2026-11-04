#!/usr/bin/env python3
# The simple cases' decks (Chunhui Zhao, 2026-10-31: "go ahead with the simple
# cases, run them on local Mac"; plan "2D plane strain simple tests plan"):
# twelve decks made from config/fd_cycle/mapview_rs_fd.toml, written to
# simple_cases/configs/<case>.toml; a thirteenth, D5f (2026-10-31: "let's
# prepare for this case, a smooth design around this bend, and rerun"), is
# D1f with its 1.5 deg corner rounded into a 10 km arc; two more (2026-10-31:
# "let's refine the mesh 1 level / 2 levels for fault and rerun this
# problem") are D5f with 250 m and 125 m on the fault; three more (P0 of
# vertex-normal-plan-rev1, 2026-11-01: "send this document to 2d thread and
# start implementation and test") are D5f's event 1 with the Nitsche penalty
# scaled by 1, 2 and 4 (OPTS below); five more (P3 and P4 of the plan,
# 2026-11-01) turn on [fault_geometry] vertex_normals = "averaged": D0's and
# D0-averaged's event 1 (the straight control off and on), D5f-averaged's
# event 1 (the smoke run), and D5f and its 250 m refinement averaged for
# 12,000 yr; five more (2026-11-01: "do you suggest we run a few more cases
# by adjusting the angles (similar to what UCERF fault)"; "we can smooth the
# trace bends if needed") are smoothed bends at 250 m with averaged normals:
# 3 and 8 deg over 10 km and over 2.5 km arcs, and a +3 / -3 deg pair of
# 5 km arcs (SMOOTH below); three more (2026-11-01: "add functionality
# about velocity strengthening patch near the tip of the fault, and test on
# previous 1.5 degree fault with averaged normal") are P4's 250 m D5f with
# [cycle.fields.a] kind = "tip_zone" (a measured along the fault from each
# tip): stronger (a = 0.08), wider (from 6 km of each tip) and both
# (TIPZONE below); two more (2026-11-02: "implement the relaxation, test and
# verify it, and apply to UCERF 2D plane strain, with 10MPa cap") are D9f
# with the floor at 10 MPa for 30,000 yr, without and with sigma_n'
# relaxation, t_r 5 kyr (RELAX below); five more (P3 of the PML design plan
# rev2, test T5, 2026-11-04: "start P3") are D0 to the end of event 1 and its
# hand-back, 900 yr: the absorbing walls (the reference), the CFS-PML band at
# f_alpha = 0.005, 0.01 and 0.02 Hz, and the sponge comparator (PML below).
# Run from miniapps/seas:
#   python3 simple_cases/scripts/make_decks.py [--check]
# --check: write nothing, exit 1 if a committed deck differs from what this
# script makes (the decks are committed as made).
#
# Edits to mapview_rs_fd.toml, every case:
#   [mesh] path        simple_cases/mesh/<mesh>.msh (+-240 km box, 500 m)
#   [time] tfinal      12,000 yr of 365 days (about 3 events, see README)
#   [cycle.output]     prefix = <case>; [[cycle.output.station]] on the
#                      fault by arc length (the trace's own points, so a
#                      bent limb's stations sit on it), plus the bend set
#   [cycle.benchmark]  region_x3 = [-1e4, 1e4] (the bent limb reaches
#                      y = 1.31 km), as the bend convergence decks
#   [cycle.snapshots.fault] dt_dyn = 0.05 s (was 0.1 s; output only)
#   [cycle.snapshots.volume] dt_dyn = 1.0 s (was 0.5 s; 1 MB per snapshot,
#                      so about 0.27 GB per fully dynamic run of 3 events)
# QD only (Q*):        [cycle] switching = false
# fully dynamic (D*):  [cycle.output] handback_table = true
# floor (D1f, D2f, D5f): [cycle] sigma_n_floor_pa = 1.0e6 (UCERF run8's)
# P0 penalty scan (OPTS): [cycle] penalty_factor = 1, 2, 4; [numerics]
#                      cfl = 0.25 in all three (D5f: 0.5); [cycle.switch]
#                      max_dynamic_phases = 1 (stop at the end of event 1)
# vertex-normal P3/P4 (OPTS): [fault_geometry] vertex_normals = "averaged"
#                      (the "_avg" cases; turn threshold 5 deg, the default);
#                      the "_ev1" cases max_dynamic_phases = 1
# smoothed bends (SMOOTH): as P4's 250 m case (averaged, floor), on their
#                      own meshes; region_x3's note gives the fault's own |y|
# tip zones (TIPZONE): as P4's 250 m case (averaged, floor, its mesh), with
#                      [cycle.fields.a] kind = "tip_zone" in place of the
#                      rect taper; region_x3's note gives the fault's own |y|
# relaxation (RELAX):  as D9f (averaged, its mesh), with [cycle]
#                      sigma_n_floor_pa = 10.0e6, [time] tfinal 30,000 yr,
#                      and in the _relax case [cycle.relaxation]
#                      t_relax_yr = 5000
# PML hand-back (PML):  as D0, [time] tfinal 900 yr (event 1 at ~827 yr, its
#                      hand-back, then QD); the band cases [cycle.absorbing]
#                      faces = "off" and [cycle.pml] on every outer face
#                      (attributes 5 and 1), 150 km (30 elements of 5 km),
#                      N 2, R0 1e-3; kind "cfs" with f_alpha_hz and
#                      wall_speed_factor = "auto" (D3), or kind "sponge";
#                      P4 (T6): the relax5k D9f deck with the band at
#                      f_alpha 0.01 Hz, 12,000 yr (its first two events)
import argparse
import math
import os
import sys

BASE = "config/fd_cycle/mapview_rs_fd.toml"
YR = 365.0 * 86400.0
T_FINAL_YR = 12000.0
FAULT_HALF = 15.0e3

# mesh: (theta_deg, arc [m], s_bend [m], lc_fault [m]), as the README's gmsh
# lines and bend.geo; lc_fault is the fault points' spacing in arc length
MESHES = {
    "straight_h500": (0.0, 0.0, 0.0, 500.0),
    "sharp015_h500": (1.5, 0.0, 0.0, 500.0),
    "sharp05_h500": (5.0, 0.0, 0.0, 500.0),
    "arc05L10_h500": (5.0, 10000.0, 0.0, 500.0),
    "sharp015vs_h500": (1.5, 0.0, 12500.0, 500.0),
    "arc015L10_h500": (1.5, 10000.0, 0.0, 500.0),
    "arc015L10_h250": (1.5, 10000.0, 0.0, 250.0),
    "arc015L10_h125": (1.5, 10000.0, 0.0, 125.0),
    # the smoothed bends (2026-11-01), 250 m: 3 and 8 deg over 10 and 2.5 km
    "arc03L10_h250": (3.0, 10000.0, 0.0, 250.0),
    "arc08L10_h250": (8.0, 10000.0, 0.0, 250.0),
    "arc03L2p5_h250": (3.0, 2500.0, 0.0, 250.0),
    "arc08L2p5_h250": (8.0, 2500.0, 0.0, 250.0),
    "sbend03L05_h250": (3.0, 5000.0, -3750.0, 250.0),
}
# a second, smooth bend (bend.geo's theta2_deg, arc2, s_bend2): +3 deg over
# 5 km centred at s = -3.75 km, then -3 deg over 5 km centred at +3.75 km
BEND2 = {"sbend03L05_h250": (-3.0, 5000.0, 3750.0)}

# case: (mesh, regime, floor, plan row), in the run order
CASES = [
    ("D1f_sharp015_floor", "sharp015_h500", "D", True, "D1 with the 1 MPa floor: a 1.5 deg sharp bend at x = 0"),
    ("D2f_sharp05_floor", "sharp05_h500", "D", True, "D2 with the 1 MPa floor: a 5 deg sharp bend at x = 0"),
    ("Q1_sharp015", "sharp015_h500", "Q", False, "Andy's [2]: UCERF's angle, a 1.5 deg sharp bend at x = 0"),
    ("Q0_straight", "straight_h500", "Q", False, "Andy's [1]: the straight control"),
    ("Q2_sharp05", "sharp05_h500", "Q", False, "Andy's [2]: a 5 deg sharp bend at x = 0"),
    ("Q3_arc05L10", "arc05L10_h500", "Q", False, "Andy's [2]: 5 deg over a 10 km arc centred at x = 0"),
    ("Q4_sharp015vs", "sharp015vs_h500", "Q", False, "addition: a 1.5 deg sharp bend at x = 12.5 km, in the creeping part"),
    ("D0_straight", "straight_h500", "D", False, "Andy's [3]: the straight control"),
    ("D1_sharp015", "sharp015_h500", "D", False, "Andy's [4]: UCERF's angle, a 1.5 deg sharp bend at x = 0"),
    ("D2_sharp05", "sharp05_h500", "D", False, "Andy's [4]: a 5 deg sharp bend at x = 0"),
    ("D3_arc05L10", "arc05L10_h500", "D", False, "Andy's [4]: 5 deg over a 10 km arc centred at x = 0"),
    ("D4_sharp015vs", "sharp015vs_h500", "D", False, "addition: a 1.5 deg sharp bend at x = 12.5 km, in the creeping part"),
    ("D5f_arc015L10_floor", "arc015L10_h500", "D", True, "D1f with the bend rounded: 1.5 deg over a 10 km arc centred at x = 0"),
    ("D5f_arc015L10_h250_floor", "arc015L10_h250", "D", True, "D5f with the fault mesh refined once: 250 m on the fault"),
    ("D5f_arc015L10_h125_floor", "arc015L10_h125", "D", True, "D5f with the fault mesh refined twice: 125 m on the fault"),
    ("D5f_arc015L10_pf1_ev1", "arc015L10_h500", "D", True, "P0 of the vertex-normal plan: D5f's event 1, penalty_factor 1, cfl 0.25"),
    ("D5f_arc015L10_pf2_ev1", "arc015L10_h500", "D", True, "P0 of the vertex-normal plan: D5f's event 1, penalty_factor 2, cfl 0.25"),
    ("D5f_arc015L10_pf4_ev1", "arc015L10_h500", "D", True, "P0 of the vertex-normal plan: D5f's event 1, penalty_factor 4, cfl 0.25"),
    ("D0_straight_ev1", "straight_h500", "D", False, "P3 of the vertex-normal plan: D0's event 1, face normals (the regression's off run)"),
    ("D0_straight_avg_ev1", "straight_h500", "D", False, "P3 of the vertex-normal plan: D0's event 1 with averaged vertex normals (on a straight fault: as off)"),
    ("D5f_arc015L10_avg_ev1", "arc015L10_h500", "D", True, "P3 of the vertex-normal plan: D5f's event 1 with averaged vertex normals (the smoke run)"),
    ("D5f_arc015L10_avg_floor", "arc015L10_h500", "D", True, "P4 of the vertex-normal plan: D5f with averaged vertex normals"),
    ("D5f_arc015L10_h250_avg_floor", "arc015L10_h250", "D", True, "P4 of the vertex-normal plan: D5f at 250 m with averaged vertex normals"),
    # the smoothed bends, in the run order (8 deg first)
    ("D7f_arc08L10_h250_avg_floor", "arc08L10_h250", "D", True, "smoothed bend: 8 deg over a 10 km arc (R 71.6 km), averaged vertex normals"),
    ("D9f_arc08L2p5_h250_avg_floor", "arc08L2p5_h250", "D", True, "smoothed bend: 8 deg over a 2.5 km arc (R 17.9 km), averaged vertex normals"),
    ("D6f_arc03L10_h250_avg_floor", "arc03L10_h250", "D", True, "smoothed bend: 3 deg over a 10 km arc (R 191 km), averaged vertex normals"),
    ("D8f_arc03L2p5_h250_avg_floor", "arc03L2p5_h250", "D", True, "smoothed bend: 3 deg over a 2.5 km arc (R 47.7 km), averaged vertex normals"),
    ("D10f_sbend03L05_h250_avg_floor", "sbend03L05_h250", "D", True, "smoothed bends: +3 deg then -3 deg, 5 km arcs (R 95.5 km) 2.5 km apart, averaged vertex normals"),
    # the tip zones on P4's 250 m D5f (2026-11-01)
    ("D5f_arc015L10_h250_avg_floor_T1a08", "arc015L10_h250", "D", True, "tip zone T1 (stronger): P4's 250 m D5f with a = 0.08 within 3 km of each tip, 2 km taper"),
    ("D5f_arc015L10_h250_avg_floor_T2w8", "arc015L10_h250", "D", True, "tip zone T2 (wider): P4's 250 m D5f with a = 0.04 within 6 km of each tip, 2 km taper"),
    ("D5f_arc015L10_h250_avg_floor_T3a08w8", "arc015L10_h250", "D", True, "tip zone T3 (both): P4's 250 m D5f with a = 0.08 within 6 km of each tip, 2 km taper"),
    ("D9f_arc08L2p5_h250_avg_floor10", "arc08L2p5_h250", "D", True, "relaxation control: D9f with the sigma_n' floor at 10 MPa, 30,000 yr, no relaxation"),
    ("D9f_arc08L2p5_h250_avg_floor10_relax5k", "arc08L2p5_h250", "D", True, "relaxation V1: D9f with sigma_n' relaxation (t_r 5 kyr) and the floor at 10 MPa, 30,000 yr"),
    # PML plan rev2 P3 (test T5): event 1 and its hand-back, 900 yr
    ("D0_straight_hb_ref", "straight_h500", "D", False, "PML P3 (T5) reference: D0's event 1 and hand-back with the absorbing walls"),
    ("D0_straight_hb_pml0p005", "straight_h500", "D", False, "PML P3 (T5): D0's event 1 and hand-back, CFS-PML band, f_alpha 0.005 Hz"),
    ("D0_straight_hb_pml0p01", "straight_h500", "D", False, "PML P3 (T5): D0's event 1 and hand-back, CFS-PML band, f_alpha 0.01 Hz"),
    ("D0_straight_hb_pml0p02", "straight_h500", "D", False, "PML P3 (T5): D0's event 1 and hand-back, CFS-PML band, f_alpha 0.02 Hz"),
    ("D0_straight_hb_sponge", "straight_h500", "D", False, "PML P3 (T5) comparator: D0's event 1 and hand-back, sponge band (same profile)"),
    # PML plan rev2 P4 (T6): the relaxation V1 deck with the CFS-PML band, 12,000 yr
    ("D9f_arc08L2p5_h250_avg_floor10_relax5k_pml0p01", "arc08L2p5_h250", "D", True, "PML P4 (T6): relaxation V1 (D9f, t_r 5 kyr, 10 MPa floor) with the CFS-PML band, f_alpha 0.01 Hz, 12,000 yr"),
]
SMOOTH = ["D7f_arc08L10_h250_avg_floor", "D9f_arc08L2p5_h250_avg_floor", "D6f_arc03L10_h250_avg_floor",
          "D8f_arc03L2p5_h250_avg_floor", "D10f_sbend03L05_h250_avg_floor"]

# P0 of vertex-normal-plan-rev1 (/mnt/project-files/2d-inplane/
# nonplanar-literature/): if the arc's vertex jump in sigma_n' is the Nitsche
# penalty acting on the two faces' slip-direction mismatch, it grows with
# penalty_factor.  The explicit step does not see the penalty
# (solver/nitsche_cfl.hpp: dt = cfl C(p) h / (p^2 c_p)); the face stiffness
# grows with it, so omega_max can grow up to about sqrt(4) = 2 at factor 4.
# cfl 0.25 keeps D5f's margin (cfl 0.5 at factor 1) at factor 4; the factor
# 1 run at the same cfl separates the step's own effect from the penalty's.
OPTS = {
    "D5f_arc015L10_pf1_ev1": dict(penalty=1.0, cfl=0.25, phases=1),
    "D5f_arc015L10_pf2_ev1": dict(penalty=2.0, cfl=0.25, phases=1),
    "D5f_arc015L10_pf4_ev1": dict(penalty=4.0, cfl=0.25, phases=1),
    # vertex-normal plan P3 (regression and smoke) and P4 (the D5f runs)
    "D0_straight_ev1": dict(phases=1),
    "D0_straight_avg_ev1": dict(phases=1, averaged=True),
    "D5f_arc015L10_avg_ev1": dict(phases=1, averaged=True),
    "D5f_arc015L10_avg_floor": dict(averaged=True),
    "D5f_arc015L10_h250_avg_floor": dict(averaged=True),
}
for _c in SMOOTH:
    OPTS[_c] = dict(averaged=True)

# the tip zones: [cycle.fields.a] kind = "tip_zone" (cycle/fault_chains.hpp),
# d the distance along the fault to the nearer free tip: a = tip for
# d <= width - taper, linear to inside at d = width, inside beyond.  The
# rect taper of mapview_rs_fd.toml (0.004 for |x| <= 10 km, 0.04 from 12 km)
# is width 5 km, taper 2 km, tip 0.04 on this 30 km fault (to the bent
# limb's x, which is 5 m shorter than its arc length at the east tip).
# b = 0.03 everywhere; the prestress takes the local a ([cycle.prestress] a
# = "local"), so the zone starts at steady state at V_p for its own a.
TIPZONE = {
    "D5f_arc015L10_h250_avg_floor_T1a08": dict(width=5.0e3, taper=2.0e3, tip=0.08),
    "D5f_arc015L10_h250_avg_floor_T2w8": dict(width=8.0e3, taper=2.0e3, tip=0.04),
    "D5f_arc015L10_h250_avg_floor_T3a08w8": dict(width=8.0e3, taper=2.0e3, tip=0.08),
}
for _c in TIPZONE:
    OPTS[_c] = dict(averaged=True)

# sigma_n' relaxation (relaxation design rev 0, V1): D9f, the tightest
# smoothed bend (R 17.9 km, 25 MPa sigma_n'), whose two events in 12,000 yr
# put 3 DOFs on the 1 MPa floor at event 2; here the floor is 10 MPa (the
# UCERF rerun's) and the run 30,000 yr (about 5 events), once without and
# once with relaxation ([cycle.relaxation], cycle/normal_stress_relaxation.hpp).
RELAX = {
    "D9f_arc08L2p5_h250_avg_floor10": dict(floor=10.0e6, years=30000.0, t_relax_yr=None),
    "D9f_arc08L2p5_h250_avg_floor10_relax5k": dict(floor=10.0e6, years=30000.0, t_relax_yr=5000.0),
    # PML P4: the relax5k deck with the band, its first two events
    "D9f_arc08L2p5_h250_avg_floor10_relax5k_pml0p01": dict(floor=10.0e6, years=12000.0, t_relax_yr=5000.0),
}
for _c in RELAX:
    OPTS[_c] = dict(averaged=True)


def trace(u, theta_deg, arc, s_bend, bend2=None):
    """Point (x, y) [m] of the fault at u, its arc length from the centre of
    the straight west part (u = x there), as mesh/bend.geo builds it;
    bend2 = (theta2_deg, arc2, s_bend2), the second smooth bend, or None."""
    if bend2 is not None:
        th2, arc2, sb2 = math.radians(bend2[0]), bend2[1], bend2[2]
        a2s, a2e = sb2 - arc2 / 2.0, sb2 + arc2 / 2.0
        if u > a2s:
            th, k2 = math.radians(theta_deg), th2 / arc2
            X2, Y2 = trace(a2s, theta_deg, arc, s_bend)
            if u <= a2e:
                ph = th + (u - a2s) * k2
                return X2 + (math.sin(ph) - math.sin(th)) / k2, Y2 - (math.cos(ph) - math.cos(th)) / k2
            X3 = X2 + (math.sin(th + th2) - math.sin(th)) / k2
            Y3 = Y2 - (math.cos(th + th2) - math.cos(th)) / k2
            return X3 + (u - a2e) * math.cos(th + th2), Y3 + (u - a2e) * math.sin(th + th2)
    th = math.radians(theta_deg)
    s = u - s_bend
    if arc == 0.0:
        x, y = (s, 0.0) if s <= 0.0 else (s * math.cos(th), s * math.sin(th))
    else:
        R, s0 = arc / th, -arc / 2.0
        if s <= s0:
            x, y = s, 0.0
        elif s <= -s0:
            p = (s - s0) / R
            x, y = s0 + R * math.sin(p), R * (1.0 - math.cos(p))
        else:
            x = s0 + R * math.sin(th) + (s + s0) * math.cos(th)
            y = R * (1.0 - math.cos(th)) + (s + s0) * math.sin(th)
    return x + s_bend, y


def stations(mesh):
    """Fault stations on the trace: fltst_s<km> every 3 km of arc length,
    and the bend set: bend_w and bend_e 1 m from the vertex (or the arc's
    mid point) on the west and the east face, so each face keeps its own
    sigma_n'; bend_s<+-km> along the fault from the bend; the arc's ends."""
    th, arc, sb, lc = MESHES[mesh]
    b2 = BEND2.get(mesh)
    out = []
    for k in range(-12, 13, 3):
        out.append((f"fltst_s{k:+03d}", trace(k * 1e3, th, arc, sb, b2)))
    # the bend set, and with a second bend the same set again (bend2_...)
    for tag, c, a in [("bend", sb, arc)] + ([("bend2", b2[2], b2[1])] if b2 else []):
        # 1 m along the face's chord (the fault is straight segments of lc
        # between points on the curve), so the point lies on that face
        v = trace(c, th, arc, sb, b2)
        for name, u in ((tag + "_w", c - lc), (tag + "_e", c + lc)):
            n = trace(u, th, arc, sb, b2)
            L = math.hypot(n[0] - v[0], n[1] - v[1])
            out.append((name, (v[0] + (n[0] - v[0]) / L, v[1] + (n[1] - v[1]) / L)))
        offs = [-3.0, -1.0, 1.0, 3.0]
        if c + 3e3 > FAULT_HALF:
            offs = [-3.0, -1.0, 1.0, 2.0]       # the tip is 2.5 km past the bend
        if a > 0.0:
            offs = [-a / 2e3] + offs + [a / 2e3]
        for d in offs:
            # +-1.25 km (a 2.5 km arc's ends) keeps its second decimal
            nm = f"{d:+.1f}" if abs(d * 10.0 - round(d * 10.0)) < 1e-9 else f"{d:+.2f}"
            out.append((f"{tag}_s{nm}", trace(c + d * 1e3, th, arc, sb, b2)))
    return out


# the request each deck answers (Chunhui Zhao, 2026-10-31)
GO = "# 2026-10-31: \"go ahead with the simple cases, run them on local Mac\"):"
GO_CASE = {"D5f_arc015L10_floor":
           "# 2026-10-31: \"a smooth design around this bend, and rerun\"):",
           "D5f_arc015L10_h250_floor":
           "# 2026-10-31: \"let's refine the mesh 1 level / 2 levels for fault and rerun\"):",
           "D5f_arc015L10_h125_floor":
           "# 2026-10-31: \"let's refine the mesh 1 level / 2 levels for fault and rerun\"):"}
for _c in OPTS:
    GO_CASE[_c] = "# 2026-11-01: \"send this document to 2d thread and start implementation and test\"):"
for _c in SMOOTH:
    GO_CASE[_c] = "# 2026-11-01: \"we can smooth the trace bends if needed\"):"
for _c in TIPZONE:
    GO_CASE[_c] = "# 2026-11-01: \"add functionality about velocity strengthening patch\"):"
for _c in RELAX:
    GO_CASE[_c] = "# 2026-11-02: \"implement the relaxation, test and verify it\"):"

# PML plan rev2 P3 (T5): kind None = the absorbing walls (the reference)
PML = {
    "D0_straight_hb_ref": dict(kind=None),
    "D0_straight_hb_pml0p005": dict(kind="cfs", f_alpha=0.005),
    "D0_straight_hb_pml0p01": dict(kind="cfs", f_alpha=0.01),
    "D0_straight_hb_pml0p02": dict(kind="cfs", f_alpha=0.02),
    "D0_straight_hb_sponge": dict(kind="sponge"),
    # P4 (T6): against D9f_arc08L2p5_h250_avg_floor10_relax5k's first two events
    "D9f_arc08L2p5_h250_avg_floor10_relax5k_pml0p01": dict(kind="cfs", f_alpha=0.01, years=12000.0),
}
PML_YEARS = 900.0
PML_L_KM = 150.0
for _c in PML:
    GO_CASE[_c] = "# 2026-11-04: \"start P3\" of the PML design plan rev2, test T5):"
GO_CASE["D9f_arc08L2p5_h250_avg_floor10_relax5k_pml0p01"] = (
    "# 2026-11-04: \"yes let's test D9f first, then UCERF\", PML plan rev2 P4):")


def header(case, mesh, regime, floor, row):
    lines = [
        f"# simple_cases/configs/{case}.toml — the simple cases (Chunhui Zhao,",
        GO_CASE.get(case, GO),
        f"# {row}.",
        f"# {'QD only' if regime == 'Q' else 'Fully dynamic (QD and dynamic regimes)'}, mesh {mesh}"
        + ((f", sigma_n' floor {RELAX[case]['floor'] / 1e6:g} MPa." if case in RELAX
            else ", sigma_n' floor 1 MPa.") if floor else ", no sigma_n' floor."),
        "# Made by simple_cases/scripts/make_decks.py from",
        "# config/fd_cycle/mapview_rs_fd.toml; the edits are listed there and in",
        "# simple_cases/README.md.  Plane strain, no 2.5D term, +-240 km box, all",
        (f"# walls Dirichlet with a {PML_L_KM:.0f} km {'sponge' if PML[case]['kind'] == 'sponge' else 'CFS-PML'} band, "
         + (f"{PML[case]['years']:,.0f} yr." if "years" in PML[case] else f"{PML_YEARS:.0f} yr (event 1, its hand-back).")
         if PML.get(case, {}).get("kind")
         else f"# walls absorbing in the dynamic regime, {PML_YEARS:.0f} yr (event 1, its hand-back)."
         if case in PML
         else "# walls absorbing in the dynamic regime, to the end of event 1."
         if OPTS.get(case, {}).get("phases") == 1
         else f"# walls absorbing in the dynamic regime, {RELAX.get(case, {}).get('years', T_FINAL_YR):,.0f} yr."),
        "#",
        "# Run (Mac, from miniapps/seas): simple_cases/scripts/run_simple_mac.sh",
        "#",
        "# mapview_rs_fd.toml's notes follow (its mesh: the straight fault, +-60 km).",
        "#",
    ]
    return "\n".join(lines) + "\n"


def make(case, mesh, regime, floor, row, base):
    t = base

    def sub(old, new, count=1):
        nonlocal t
        assert t.count(old) == count, (case, old, t.count(old))
        t = t.replace(old, new)

    sub('path  = "mapview_rs/mesh/mapview_rs_500m.msh"', f'path  = "simple_cases/mesh/{mesh}.msh"')
    lc = MESHES[mesh][3]
    if lc != 500.0:
        sub("# 500 m on the fault, in m.", f"# {lc:.0f} m on the fault (mapview_rs_fd.toml: 500 m), in m.")
        sub("# Dense G: N_f = 120 at p = 1 on the 500 m mesh.",
            f"# Dense G: N_f = {round(2 * 30e3 / lc)} at p = 1 on this {lc:.0f} m mesh (120 at 500 m).")
    years = (PML[case].get("years", PML_YEARS) if case in PML
             else RELAX.get(case, {}).get("years", T_FINAL_YR))
    sub("# 1500 yr of 365 days.\n[time]\ntfinal     = 4.7304e10",
        f"# {years:,.0f} yr of 365 days (mapview_rs_fd.toml: 1500 yr).\n[time]\n"
        f"tfinal     = {years * YR:.6e}")
    if regime == "Q":
        sub("switching            = true", "switching            = false   # QD only (mapview_rs_fd.toml: true)")
    if floor and case in RELAX:
        sub("penalty_factor       = 1.0       # both regimes' K\n",
            "penalty_factor       = 1.0       # both regimes' K\n"
            f"sigma_n_floor_pa     = {RELAX[case]['floor'] / 1e6:.1f}e6    # the UCERF rerun's 10 MPa floor (Chunhui Zhao, 2026-11-02:"
            " \"with 10MPa cap\"; D9f: 1 MPa)\n")
    elif floor:
        sub("penalty_factor       = 1.0       # both regimes' K\n",
            "penalty_factor       = 1.0       # both regimes' K\n"
            "sigma_n_floor_pa     = 1.0e6     # the 1 MPa floor of UCERF run8 (not in mapview_rs_fd.toml)\n")
    if RELAX.get(case, {}).get("t_relax_yr"):
        sub("# The BP5-FD switch.\n[cycle.switch]\n",
            "# sigma_n' relaxation (make_decks.py RELAX; cycle/normal_stress_relaxation.hpp):\n"
            "# the friction law sees sigma0 + dsig - r, dr/dt = (dsig - r) / t_r in the QD\n"
            "# regime, r frozen in the dynamic regime; the outputs keep the elastic sigma_n'.\n"
            "[cycle.relaxation]\n"
            f"t_relax_yr = {RELAX[case]['t_relax_yr']:.1f}\n\n"
            "# The BP5-FD switch.\n[cycle.switch]\n")
    pk = PML.get(case, {}).get("kind")
    if pk:
        sub('[cycle.absorbing]\nfaces = "all"\n',
            '[cycle.absorbing]\nfaces = "off"   # the PML band instead (mapview_rs_fd.toml: "all")\n\n'
            "# PML design plan rev2 (cycle/cycle_config.hpp PmlConfig): the band on every\n"
            "# outer face, the walls Dirichlet in both regimes"
            + (".  CFS-PML: QD with the static stretch\n"
               "# (K0), the dynamic phase in the time domain; the wall speed times\n"
               "# (W + 2 Delta L) / W (D3).\n" if pk == "cfs" else
               ".  Sponge (the comparator):\n"
               "# the dynamic damping rho sum_a d_a U' alone, the same profile; QD and the\n"
               "# walls unchanged.\n")
            + "[cycle.pml]\n"
            "on           = true\n"
            f'kind         = "{pk}"\n'
            "wall_attrs   = [5, 1]\n"
            f"thickness_km = {PML_L_KM:.1f}   # 30 elements of 5 km\n"
            "power        = 2\n"
            "R0           = 1.0e-3\n"
            + (f"f_alpha_hz   = {PML[case]['f_alpha']}\n"
               'wall_speed_factor = "auto"\n' if pk == "cfs" else ""))
    o = OPTS.get(case)
    if o:
        if o.get("penalty", 1.0) != 1.0:
            sub("penalty_factor       = 1.0       # both regimes' K\n",
                f"penalty_factor       = {o['penalty']:.1f}       # both regimes' K (mapview_rs_fd.toml: 1.0)\n")
        if "cfl" in o:
            sub("cfl       = 0.5\n", f"cfl       = {o['cfl']}      # P0: penalty up to x 4 (mapview_rs_fd.toml: 0.5)\n")
        if "phases" in o:
            why = "P0" if "penalty" in o else "P3"
            sub("max_dynamic_phases = 0\n",
                f"max_dynamic_phases = {o['phases']}   # {why}: stop at the end of event 1 (mapview_rs_fd.toml: 0)\n")
        if o.get("averaged"):
            sub("pin_min_margin  = 0.1\n",
                "pin_min_margin  = 0.1\n"
                "vertex_normals  = \"averaged\"   # "
                + ("averaged vertex normals" if case in SMOOTH or case in TIPZONE or case in RELAX
                   else "vertex-normal plan P3/P4")
                + " (mapview_rs_fd.toml: \"face\")\n")
    sub('prefix          = "mapview_rs_fd"', f'prefix          = "{case}"'
        + ("\nhandback_table  = true" if regime == "D" else ""))
    # stations: replace the nine fltst_x entries
    i0 = t.index("# On-fault stations every 3 km from -12 to 12 km.")
    i1 = t.index("# Off-fault stations")
    st = ["# On-fault stations on the trace (make_decks.py): fltst_s<km> every 3 km",
          "# of arc length from -12 to 12 km; bend_w / bend_e 1 m from the bend",
          "# vertex (arc: its mid point) on the west / east face; bend_s<km> along",
          "# the fault from the bend (arc: also its two ends)."]
    for name, (x, y) in stations(mesh):
        st += ["[[cycle.output.station]]", f'name = "{name}"', f"x    = [{x!r}, {y!r}]"]
    t = t[:i0] + "\n".join(st) + "\n\n" + t[i1:]
    if case in TIPZONE:
        z = TIPZONE[case]
        sub("# a: 0.004 for |x| <= 10 km, rising linearly to 0.04 at 12 km.\n"
            "[cycle.fields.a]\n"
            "center  = [0.0, 0.0]\n"
            "axis1   = [1.0, 0.0]\n"
            "half1   = 10.0e3\n"
            "taper   = 2.0e3\n"
            "inside  = 0.004\n"
            "outside = 0.04\n",
            f"# a (make_decks.py TIPZONE; mapview_rs_fd.toml: the rect taper, 0.004 for\n"
            f"# |x| <= 10 km, rising linearly to 0.04 at 12 km): measured along the fault\n"
            f"# from each free tip, d: a = {z['tip']} for d <= {(z['width'] - z['taper']) / 1e3:g} km, "
            f"linear to 0.004 at\n"
            f"# d = {z['width'] / 1e3:g} km, 0.004 beyond (b = 0.03).\n"
            "[cycle.fields.a]\n"
            "kind   = \"tip_zone\"\n"
            f"width  = {z['width'] / 1e3:.1f}e3\n"
            f"taper  = {z['taper'] / 1e3:.1f}e3\n"
            "inside = 0.004\n"
            f"tip    = {z['tip']}\n")
    if case in SMOOTH or case in TIPZONE or case in RELAX:
        th, arc, sb, _ = MESHES[mesh]
        ymax = max(abs(trace(k * 250.0, th, arc, sb, BEND2.get(mesh))[1]) for k in range(-60, 61))
        sub("region_x3       = [-1.0, 1.0]", f"region_x3       = [-1.0e4, 1.0e4]   # the fault reaches |y| = {ymax / 1e3:.2f} km")
    else:
        sub("region_x3       = [-1.0, 1.0]", "region_x3       = [-1.0e4, 1.0e4]   # the bent limb reaches y = 1.31 km")
    sub("dt_dyn = 0.1\n", "dt_dyn = 0.05           # mapview_rs_fd.toml: 0.1 (sigma_n' against slip at the bend)\n")
    sub("dt_dyn  = 0.5\n", "dt_dyn  = 1.0           # mapview_rs_fd.toml: 0.5 (1 MB each; disk)\n")
    return header(case, mesh, regime, floor, row) + t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    base = open(BASE).read()
    bad = 0
    for case, mesh, regime, floor, row in CASES:
        text = make(case, mesh, regime, floor, row, base)
        path = f"simple_cases/configs/{case}.toml"
        if a.check:
            same = os.path.exists(path) and open(path).read() == text
            print(f"{'same' if same else 'DIFFERS'} {path}")
            bad += not same
        else:
            with open(path, "w") as f:
                f.write(text)
            print("wrote", path)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
