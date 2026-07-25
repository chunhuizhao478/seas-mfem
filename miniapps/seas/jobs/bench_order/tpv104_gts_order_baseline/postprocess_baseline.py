#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TPV104 200 m GTS order baseline — log/Caliper post-processor.

Parses the SIX legs of the speed baseline

    MFEM  p1 / p2 / p3   (seas_spatial_dyn_driver, ADER order p+1, [numerics].lts="off")
    SeisSol o2 / o3 / o4 (compile-time ConvergenceOrder, ClusteredLTS = 1)

on the SCEC TPV104 200 m production mesh (2,464,689 tets, md5-identical for both
codes), t = 2.0 s, Courant 0.5, and emits the baseline table:

  * per-leg  seconds of wall per simulated second   (s / sim-s)
  * per-leg  simulated seconds per wall hour        (sim-s / wall-h)
  * per-leg  per-element-update cost in MICROSECONDS x CORE
  * per-leg  steps, dt, elements, cores
  * the MFEM / SeisSol ratio at each ACCURACY-MATCHED order  p_k <-> o_(k+1)
  * the within-code p1 / p2 / p3 scaling, with the exactly-known (2p+1) CFL
    step-count penalty divided out so per-step work is legible on its own
  * the fairness gates, and an explicit list of everything it could NOT compute
  * a one-line CAVEAT of what this baseline cannot claim

Standard library only (no numpy / pandas / toml).  Python 3.6+.

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------
    python3 postprocess_baseline.py <RUN_DIR> [options]

<RUN_DIR> is scanned RECURSIVELY.  Leg identification is CONTENT-FIRST (the
logs identify themselves), so any sane directory layout works; the conventional
one is

    <RUN_DIR>/mfem_p1/run_mfem_p1.log
    <RUN_DIR>/mfem_p1/run_mfem_p1.cali-region-report.txt
    <RUN_DIR>/mfem_p2/...            <RUN_DIR>/mfem_p3/...
    <RUN_DIR>/seissol_o2/tpv104-ss-o2.<jobid>.out
    <RUN_DIR>/seissol_o3/...         <RUN_DIR>/seissol_o4/...

Options:
    --tfinal S            simulated window in seconds        (default 2.0)
    --elements N          mesh element count                 (default 2464689)
    --dt-ref-p3 S         pinned p3/o4 dt; p1,p2 derived x7/3, x7/5
                                                             (default 3.66128e-4)
    --mfem-region NAME    Caliper region for the MFEM time loop
                                      (default seas::spatial_dyn::step)
    --mfem-cores N        override MFEM allocated cores per leg
    --seissol-cores N     override SeisSol ALLOCATED cores per leg (default 256)
    --seissol-cpus-per-rank N   cores allocated per SeisSol rank  (default 16)
    --dt-tol REL          relative dt agreement tolerance      (default 1e-5)
    --json PATH           also write the parsed/derived numbers as JSON
    --strict              treat every FAIL gate as fatal (rc 3)

Exit codes:  0 = all six legs complete and every gate PASSED
             1 = a leg (or a required quantity) is missing / uncomputable
             2 = usage error
             3 = all legs present but a fairness gate FAILED

--------------------------------------------------------------------------
WHERE EVERY NUMBER COMES FROM  (grounded; do not "improve" without re-reading)
--------------------------------------------------------------------------
MFEM  (miniapps/seas)
  T_loop   Caliper region `seas::spatial_dyn::step`, MAX-across-ranks column of
           the runtime-report.  The scope opens at
           drivers/spatial_dyn_driver.cpp:4874 as the first statement of the GTS
           step-loop body.  MFEM_PERF_SCOPE is a no-op without libcaliper
           (../../general/annotation.hpp:24 vs :31) and the driver has NO other
           wall-clock instrument -> a non-Caliper binary yields NO T_loop, and
           the SLURM elapsed is NOT a substitute (it carries a 132 MB .msh read
           + serial partition + setup, which biases the cheap p1 leg most).
  steps    Caliper region.count ("Calls") cross-checked against the last
           `step k/N` progress print (spatial_dyn_driver.cpp:5227) and the
           planned `[time] nsteps` (:3142).
  dt       `[time] dt     = <x> s` (spatial_dyn_driver.cpp:3140-3142).
  order    `fe order:` / `ader order:` banner (:1429, :1442).
  ranks    `ranks:` banner (:1465).  Pure MPI, --cpus-per-task=1 -> cores=ranks.
  safety   `cfl_dg_safety = <x>` on the `[time] cfl = ...` line (:3150-3153).
           The deck DEFAULT is 3.0 (spatial/code/spatial_friction.hpp:150); a leg
           at 3.0 takes 3x the steps and understates MFEM by 3x.  Hard gate.
  IO       `ParaView output: OFF` banner line (:3958-3961).  ParaView writes,
           per-step SCEC station writes and the in-loop checkpoint all live
           INSIDE the timed region (:5121, :5126, :5209) and CANNOT be
           subtracted afterwards -> R4 is satisfied by IO being OFF, not by
           subtraction.

SeisSol  (/Users/chunhuizhao/projects/SeisSol, v1.3.1-2135 == the Expanse binary)
  T_loop   terminal `Simulation time (compute):` line, per-rank MEAN
           (src/Solver/Simulator.cpp:150; accumulated via Stopwatch::pause,
           printed by Stopwatch::print, src/Monitoring/Stopwatch.cpp:72-79).
           The PER-EPOCH `Time spent this epoch (compute):` line is WRONG for
           this purpose: Simulator.cpp:130 prints `computeStopwatch.split()`,
           and split() = now - startTime_ (Stopwatch.cpp:35-40), so it spans the
           intervening ioStopwatch block.  This parser never reads it.
  dt       `Minimum timestep: <x> <prefix>s`
           (src/Initializer/TimeStepping/ClusterLayout.cpp:52-55, SI-prefixed by
           SIUnit::formatPrefix, src/Monitoring/Unit.cpp).
  steps    NOT printed as a total; derived as ceil(tfinal/dt) and cross-checked
           against the last `Max cluster / LTS cycle updates since sync:` line
           (src/Solver/TimeStepping/TimeCluster.cpp:887-894, every 100 steps
           under GTS).
  elements `Cell count: 2'464'689` (ClusterLayout.cpp:39; formatInteger inserts
           `'` thousands separators, Unit.cpp:33-60) and the raw integer in the
           cluster-histogram row `0 : 2464689 , 37904` (ClusterLayout.cpp:118).
  GTS      `GTS has been selected.` (LtsWeights.cpp:124) + `... 1 clustered LTS`
           in the theoretical-speedup line + a SINGLE histogram row.
  cores    `Using MPI with #ranks:` and `Using OpenMP with #threads/rank:`
           (src/SeisSol.cpp:34, :48).  The validated shape runs
           OMP_NUM_THREADS = cpus_per_task - 1 with SEISSOL_COMMTHREAD=0, so one
           core per rank is ALLOCATED BUT IDLE (16x15 = 240 busy of 256).

--------------------------------------------------------------------------
THE FAIR-METRIC FORMULAS  (recon 'fair-metric'; R3 / R6)
--------------------------------------------------------------------------
    T_sim              simulated window (2.0 s, identical both codes)      [R1]
    S                  time steps to cover T_sim = ceil(T_sim / dt)
    E                  elements (2,464,689, identical both codes)          [R1]
    C                  CORES (never ranks: 256 MFEM ranks vs 16 SeisSol
                       ranks makes any per-rank metric wrong by 16x)       [R5]
    U = S * E          element-updates in the window
    s/sim-s          = T_loop / T_sim
    sim-s/wall-h     = 3600 / (s per sim-s)
    us-core/update   = T_loop * C * 1e6 / U
    ns-core/dof-upd  = us-core/update * 1000 / b(p),  b = (p+1)(p+2)(p+3)/6
                       = 4 / 10 / 20 modes per element at p1 / p2 / p3

C is reported TWICE for every leg because the two validated rank shapes do not
normalise 1:1 (R5):
    C_alloc  = cores the job billed        (256 both codes: 2 nodes x 128)
    C_active = cores actually executing    (256 MFEM; 240 SeisSol)
The 6.7 % spread is the honest uncertainty band from the shape difference.  The
HEADLINE uses C_alloc (billing-fair, identical for both codes).

dt is NOT held equal across orders (R3): each order gets its own CFL-derived dt.
Both codes derive it from the SAME algebraic formula --
  MFEM     dt = cfl/(cfl_dg_safety*(2p+1)) * h_min / cp,  h_min = 6V/A_total
           (dynamic/wave_operator.inl:7728 + :62-119; the safety factor at
           spatial/code/spatial_friction.hpp:744-750)
  SeisSol  dt = cfl * 2*insphere / (vp * (2*ConvergenceOrder - 1))
           (src/Initializer/TimeStepping/GlobalTimestep.cpp:45-46)
-- and 2*insphere == 6V/A_total, (2O-1) == (2p+1) for O = p+1, so with
cfl_dg_safety = 1.0 the two dt are IDENTICAL at every matched order.  That is
verified here, not assumed.
"""

from __future__ import print_function

import argparse
import json
import math
import os
import re
import sys

# ---------------------------------------------------------------------------
# Pre-registered contract (recon; these are PREDICTIONS the legs must satisfy)
# ---------------------------------------------------------------------------

PINNED_ELEMENTS = 2464689          # md5-verified identical mesh, both codes
PINNED_DR_FACES = 37904
DT_REF_P3 = 3.66128e-4             # verified live, both codes, 6 digits, cfl 0.5
DEFAULT_TFINAL = 2.0
DEFAULT_REGION = "seas::spatial_dyn::step"
DEFAULT_SEISSOL_CPUS_PER_RANK = 16  # --cpus-per-task=16 in the validated sbatch

ORDERS = (1, 2, 3)                 # MFEM polynomial degree p; SeisSol o = p + 1

# Modes per tet element of an order-p modal space, x 9 elastic quantities on
# both sides.  Identical between the two codes under the p<->o+1 pairing.
MODES = {1: 4, 2: 10, 3: 20}
# Fault quadrature points per face, MFEM: FaultFaceQuadDegree = 2*order
# (dynamic/wave_operator.inl:478,665) -> MFEM triangle rules 3 / 6 / 12
# (fem/intrules.cpp:1271,1287,1304).  p3 is OVER-determined (12 QP vs 10 trace
# dofs) -- called out verbatim at drivers/spatial_dyn_driver.cpp:4211-4215.
FAULT_QP_PER_FACE = {1: 3, 2: 6, 3: 12}


def _set_refs(dt_ref, region):
    """Apply the CLI overrides to the two module-level anchors."""
    global DT_REF_P3, DEFAULT_REGION
    DT_REF_P3 = dt_ref
    DEFAULT_REGION = region


def dt_pinned(p):
    """CFL-derived dt at polynomial degree p, from the verified p3 anchor.

    dt scales EXACTLY as 1/(2p+1) on both sides, so dt(p) = dt(3) * 7/(2p+1).
    Derived (not transcribed) so the three numbers cannot drift apart.
    """
    return DT_REF_P3 * 7.0 / (2.0 * p + 1.0)


def steps_pinned(p, tfinal):
    return int(math.ceil(tfinal / dt_pinned(p)))


# ---------------------------------------------------------------------------
# SI-prefix / SeisSol formatting helpers
# ---------------------------------------------------------------------------

# src/Monitoring/Unit.cpp:22-27.  'µ' in the SeisSol source is U+00B5 (0xC2 0xB5);
# accept U+03BC and ASCII 'u' as well in case a terminal transcoded the log.
_SI_POS = {"k": 1e3, "M": 1e6, "G": 1e9, "T": 1e12, "P": 1e15,
           "E": 1e18, "Z": 1e21, "Y": 1e24, "R": 1e27, "Q": 1e30}
_SI_NEG = {"m": 1e-3, "µ": 1e-6, "μ": 1e-6, "u": 1e-6,
           "n": 1e-9, "p": 1e-12, "f": 1e-15, "a": 1e-18,
           "z": 1e-21, "y": 1e-24, "r": 1e-27, "q": 1e-30}

_NUM = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
_PFX = r"[kMGTPEZYRQ]|m|µ|μ|u|n|p|f|a|z|y|r|q"


def _si(value_str, prefix_str):
    v = float(value_str)
    if not prefix_str:
        return v
    if prefix_str in _SI_POS:
        return v * _SI_POS[prefix_str]
    if prefix_str in _SI_NEG:
        return v * _SI_NEG[prefix_str]
    return v


def parse_prefixed(text, unit="s"):
    """First `<number> [<SI prefix>]<unit>` in `text`, in base units, or None.

    Matches SIUnit::formatPrefix output, e.g. '366.1280 µs' -> 3.66128e-4,
    and the scientific fallback SIUnit::formatScientific, e.g. '3.6613e-04 s'.
    """
    m = re.search(r"(" + _NUM + r")\s*(" + _PFX + r")?" + re.escape(unit)
                  + r"(?![A-Za-z])", text)
    if not m:
        return None
    return _si(m.group(1), m.group(2) or "")


def parse_format_time(text, unit="s"):
    """Parse SIUnit::formatTime output ('3 min 42.3000 s', '1 h 5 min 2 s')."""
    total = 0.0
    hit = False
    for pat, mult in ((r"(" + _NUM + r")\s*d\b", 86400.0),
                      (r"(" + _NUM + r")\s*h\b", 3600.0),
                      (r"(" + _NUM + r")\s*min\b", 60.0)):
        m = re.search(pat, text)
        if m:
            total += float(m.group(1)) * mult
            hit = True
    # strip the d/h/min tokens before hunting the seconds term, so 'min' cannot
    # be misread and so '3 min' does not donate its 3 to the seconds parse.
    rest = re.sub(r"(" + _NUM + r")\s*(?:d|h|min)\b", " ", text)
    sec = parse_prefixed(rest, unit)
    if sec is not None:
        total += sec
        hit = True
    return total if hit else None


def parse_stopwatch(lines, label):
    """Parse one seissol::Stopwatch::print line into a dict.

    Format (src/Monitoring/Stopwatch.cpp:74-79 via utils/logger, which inserts a
    space between streamed items):

      <label> <formatTime(mean)> (per rank: (<mean> ± <std>) s ; range: [ <min> s , <max> s ])

    Returns {'mean','std','min','max','raw'} in SECONDS, or None if absent.
    The per-rank scientific triple is preferred; formatTime is the fallback.
    """
    for ln in lines:
        if label not in ln:
            continue
        out = {"raw": ln.strip()}
        m = re.search(r"\(per rank:\s*\(\s*(" + _NUM + r")\s*(?:±|\+-|\+/-)\s*("
                      + _NUM + r")\s*\)\s*s", ln)
        if m:
            out["mean"] = float(m.group(1))
            out["std"] = float(m.group(2))
        r = re.search(r"range:\s*\[\s*(" + _NUM + r")\s*s\s*,\s*(" + _NUM
                      + r")\s*s\s*\]", ln)
        if r:
            out["min"] = float(r.group(1))
            out["max"] = float(r.group(2))
        if "mean" not in out:
            head = ln.split(label, 1)[1].split("(per rank")[0]
            ft = parse_format_time(head)
            if ft is not None:
                out["mean"] = ft
                out["from_format_time"] = True
        return out if "mean" in out else None
    return None


# ---------------------------------------------------------------------------
# File discovery + classification
# ---------------------------------------------------------------------------

TEXT_EXT = (".log", ".out", ".txt", ".err", ".report")
SKIP_DIRS = {".git", "ParaView", "paraview", "__pycache__", "checkpoints"}
MAX_READ_BYTES = 96 * 1024 * 1024


def read_text(path):
    try:
        if os.path.getsize(path) > MAX_READ_BYTES:
            return None
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except (OSError, IOError):
        return None


def scan(run_dir):
    """Yield (path, text) for every plausible text artifact under run_dir."""
    for dirpath, dirnames, filenames in os.walk(run_dir):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in sorted(filenames):
            if not fn.endswith(TEXT_EXT):
                continue
            path = os.path.join(dirpath, fn)
            text = read_text(path)
            if text:
                yield path, text


def classify(path, text):
    """'mfem' | 'seissol' | 'caliper' | None.

    Order matters: an MFEM SLURM .out may embed `head -40` of a Caliper report,
    so the MFEM signature is tested first and the standalone report file is
    preferred later by the association step.
    """
    if ("seas_spatial_dyn_driver" in text
            or re.search(r"^\[time\] dt\s+=", text, re.M)
            or "fe order:" in text):
        return "mfem"
    if ("Minimum timestep" in text or "Cell count:" in text
            or "SeisSol done" in text or "GTS has been selected" in text):
        return "seissol"
    if re.search(r"^\s*Path\s", text, re.M) and "time/rank" in text.lower():
        return "caliper"
    return None


def order_from_path(path):
    """MFEM p / SeisSol o token embedded in a path, or (None, None)."""
    base = path.replace(os.sep, "/")
    mp = re.search(r"(?:^|[/_.\-])p([123])(?:[/_.\-]|$)", base)
    mo = re.search(r"(?:^|[/_.\-])o([234])(?:[/_.\-]|$)", base)
    return (int(mp.group(1)) if mp else None,
            int(mo.group(1)) if mo else None)


# ---------------------------------------------------------------------------
# MFEM parsing
# ---------------------------------------------------------------------------

def parse_mfem_log(path, text):
    d = {"path": path, "code": "mfem", "notes": []}
    lines = text.splitlines()

    def grab(pat, cast=float, flags=0):
        m = re.search(pat, text, flags | re.M)
        return cast(m.group(1)) if m else None

    d["order"] = grab(r"^fe order:\s+(\d+)", int)
    d["ader_order"] = grab(r"^ader order:\s+(\d+)", int)
    d["ranks"] = grab(r"^ranks:\s+(\d+)", int)
    d["cfl"] = grab(r"^cfl:\s+(" + _NUM + r")")
    d["tfinal_cfg"] = grab(r"^tfinal:\s+(" + _NUM + r")\s*s")
    d["mixed_flux"] = grab(r"^mixed flux:\s+(\S+)", str)
    d["time_integrator"] = grab(r"^time integrator:\s+(\S+)", str)
    d["mesh"] = grab(r"^mesh:\s+(\S+)", str)

    d["dt"] = grab(r"^\[time\] dt\s+=\s+(" + _NUM + r")\s*s")
    d["dt_cfl"] = grab(r"^\[time\] dt_cfl\s+=\s+(" + _NUM + r")\s*s")
    d["nsteps_planned"] = grab(r"^\[time\] nsteps\s+=\s+(\d+)", int)
    d["cfl_dg_safety"] = grab(r"cfl_dg_safety\s*=\s*(" + _NUM + r")")
    d["seissol_equiv_dt"] = grab(
        r"SeisSol-equivalent dt_cfl \(cfl_dg_safety=1\.0\)\s*=\s*(" + _NUM + r")\s*s")

    # Actual steps executed: the last `step k/N  t = ...  V_max = ...` print
    # (spatial_dyn_driver.cpp:5227, every 100 steps + the final step).
    last_step = None
    v_max_last = None
    v_max_peak = 0.0
    for ln in lines:
        m = re.match(r"\s*step\s+(\d+)/(\d+)\s+t\s*=\s*(" + _NUM
                     + r")\s*s\s+V_max\s*=\s*(" + _NUM + r")", ln)
        if m:
            last_step = int(m.group(1))
            d["t_last"] = float(m.group(3))
            v_max_last = float(m.group(4))
            v_max_peak = max(v_max_peak, v_max_last)
    d["steps_progress"] = (last_step + 1) if last_step is not None else None
    d["V_max_last"] = v_max_last
    d["V_max_peak"] = v_max_peak if v_max_peak > 0.0 else None

    # IO / GTS evidence
    d["paraview_off"] = ("ParaView output: OFF" in text)
    # The free-surface slice is built OUTSIDE the paraview_enabled master gate
    # (spatial_dyn_driver.cpp:3865-3871 vs the gate at :3628-3632) and is reported
    # on its OWN banner line at :3966 -- the driver happily prints "ParaView
    # output: OFF" while the slice writes every 0.05 s inside the timed region.
    # So "ParaView output: OFF" alone is NOT sufficient evidence of an IO-free wall.
    d["free_surface_banner"] = bool(re.search(r"^\s*FreeSurf:", text, re.M))
    d["lts_enabled_banner"] = ("[lts] ENABLED" in text)
    d["checkpoint_written"] = ("[checkpoint]" in text
                               and "checkpoint(s) written" in text)
    d["completed"] = bool(re.search(r"\[time\] dt\s+=", text)) and (
        last_step is not None)
    return d


# --- Caliper runtime-report -------------------------------------------------

# Column names emitted by Caliper's runtime-report controller.  Longest-first so
# the '(inc)' variants win over their exclusive-time prefixes.
_CALI_COLS = [
    "Min time/rank (inc)", "Max time/rank (inc)", "Avg time/rank (inc)",
    "Time % (inc)",
    "Min time/rank", "Max time/rank", "Avg time/rank", "Time %",
    "Calls", "Count",
]


def parse_caliper_report(path, text, region):
    """Extract the region row of a Caliper runtime-report.

    Header-driven: the header line ('Path ... Min time/rank ...') is scanned for
    known column names by POSITION, then the region's numeric tokens are mapped
    onto them positionally.  If the counts disagree the parser falls back to the
    documented default layout (Min, Max, Avg, Time %, [Calls]) and says so.

    Returns dict with keys among {'min','max','avg','pct','calls','inclusive',
    'columns','mapping_mode','row'} or None.
    """
    lines = text.splitlines()
    hdr_idx = None
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith("Path") and "time/rank" in ln:
            hdr_idx = i
            break
    header = lines[hdr_idx] if hdr_idx is not None else ""

    cols = []
    taken = []
    for name in _CALI_COLS:
        start = 0
        while True:
            pos = header.find(name, start)
            if pos < 0:
                break
            if not any(a <= pos < b for a, b in taken):
                taken.append((pos, pos + len(name)))
                cols.append((pos, name))
            start = pos + 1
    cols.sort()
    col_names = [c[1] for c in cols]

    row = None
    for ln in lines[(hdr_idx + 1) if hdr_idx is not None else 0:]:
        if region in ln:
            row = ln
            break
    if row is None:
        return None

    tail = row.split(region, 1)[1]
    toks = re.findall(_NUM, tail)
    vals = [float(t) for t in toks]

    out = {"path": path, "row": row.strip(), "columns": col_names,
           "inclusive": any("(inc)" in c for c in col_names)}

    if col_names and len(col_names) == len(vals):
        out["mapping_mode"] = "header"
        mapped = dict(zip(col_names, vals))
    else:
        out["mapping_mode"] = "positional-fallback"
        default = ["Min time/rank", "Max time/rank", "Avg time/rank", "Time %",
                   "Calls"]
        mapped = dict(zip(default, vals))

    def pick(*names):
        for n in names:
            if n in mapped:
                return mapped[n]
        return None

    out["min"] = pick("Min time/rank (inc)", "Min time/rank")
    out["max"] = pick("Max time/rank (inc)", "Max time/rank")
    out["avg"] = pick("Avg time/rank (inc)", "Avg time/rank")
    out["pct"] = pick("Time % (inc)", "Time %")
    calls = pick("Calls", "Count")
    out["calls"] = int(calls) if calls is not None else None
    return out


# ---------------------------------------------------------------------------
# SeisSol parsing
# ---------------------------------------------------------------------------

def parse_seissol_log(path, text):
    d = {"path": path, "code": "seissol", "notes": []}
    lines = text.splitlines()

    # dt: `Minimum timestep: 366.1280 µs`.  Under GTS the wiggle variants
    # ('(pre-wiggle)' / '(with wiggle ...)') never appear (ClusterLayout.cpp:58,65).
    d["dt"] = None
    for ln in lines:
        if "Minimum timestep" in ln:
            d["dt"] = parse_prefixed(ln.split("Minimum timestep", 1)[1])
            d["dt_line"] = ln.strip()
            d["dt_wiggle_variant"] = ("wiggle" in ln)
            break

    # elements: prefer the histogram row (raw int); fall back to `Cell count:`
    # whose formatInteger inserts `'` thousands separators (Unit.cpp:33-60).
    hist = []
    in_hist = False
    for ln in lines:
        if "Cluster histogram" in ln:
            in_hist = True
            continue
        if in_hist:
            m = re.search(r"(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*$", ln.rstrip())
            if m:
                hist.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
                continue
            if hist:
                in_hist = False
    d["histogram"] = hist
    d["n_clusters"] = len(hist) if hist else None
    d["elements"] = sum(h[1] for h in hist) if hist else None
    d["dr_faces"] = sum(h[2] for h in hist) if hist else None
    if d["elements"] is None:
        m = re.search(r"Cell count:\s*([0-9']+)", text)
        if m:
            d["elements"] = int(m.group(1).replace("'", ""))

    # GTS evidence
    d["gts_selected"] = ("GTS has been selected." in text)
    m = re.search(r"Theoretical speedup to GTS:\s*(" + _NUM
                  + r")\s*elementwise LTS;\s*(" + _NUM + r")\s*clustered", text)
    d["speedup_elts"] = float(m.group(1)) if m else None
    d["speedup_clts"] = float(m.group(2)) if m else None

    # exactly ONE epoch pair == no periodic output module survived (Simulator.cpp:104,112)
    d["epochs_start"] = text.count("Start simulation epoch.")
    d["epochs_end"] = text.count("End simulation epoch.")

    # timers (terminal accumulated lines ONLY -- never the per-epoch split()).
    d["t_total"] = parse_stopwatch(lines, "Simulation time (total):")
    d["t_compute"] = parse_stopwatch(lines, "Simulation time (compute):")
    d["t_io"] = parse_stopwatch(lines, "Simulation time (blocking IO):")
    d["t_init_io"] = parse_stopwatch(lines, "Time spent for initial IO:")

    # cores (SeisSol.cpp:34, :48)
    m = re.search(r"Using MPI with #ranks:\s*(\d+)", text)
    d["ranks"] = int(m.group(1)) if m else None
    m = re.search(r"Using OpenMP with #threads/rank:\s*(\d+)", text)
    d["omp_threads"] = int(m.group(1)) if m else (
        1 if "OpenMP disabled" in text else None)
    # the run sbatch echoes `=== nodes=2 ranks=16 cpus/rank=16`
    m = re.search(r"nodes=(\d+)\s+ranks=(\d+)\s+cpus/rank=(\d+)", text)
    if m:
        d["slurm_nodes"] = int(m.group(1))
        d["slurm_ranks"] = int(m.group(2))
        d["slurm_cpus_per_rank"] = int(m.group(3))

    # order: the run sbatch echoes `=== exe: .../seissol-elastic-o4-f64`
    m = re.search(r"seissol-\w+-o(\d)-f\d+", text)
    d["order"] = int(m.group(1)) if m else None

    # step audit trail (TimeCluster.cpp:892; every 100 steps under GTS)
    last = None
    for ln in lines:
        m = re.search(r"Max cluster / LTS cycle updates since sync:\s*(\d+)"
                      r"\s*at time\s*(" + _NUM + r")", ln)
        if m:
            last = (int(m.group(1)), float(m.group(2)))
    d["last_progress"] = last

    # load imbalance + FLOP (LoopStatistics.cpp:172; FlopCounter.cpp:148-168)
    m = re.search(r"Load imbalance:\s*(" + _NUM + r")\s*%", text)
    d["load_imbalance_pct"] = float(m.group(1)) if m else None
    m = re.search(r"Total calculated HW-FLOP:\s*(" + _NUM + r")\s*("
                  + _PFX + r")?FLOP", text)
    d["hw_flop"] = _si(m.group(1), m.group(2) or "") if m else None
    m = re.search(r"Total calculated NZ-FLOP:\s*(" + _NUM + r")\s*("
                  + _PFX + r")?FLOP", text)
    d["nz_flop"] = _si(m.group(1), m.group(2) or "") if m else None

    # per-element regression slopes (LoopStatistics.cpp:226-236), RANK-seconds
    d["per_element"] = {}
    for m in re.finditer(r"(compute\w+)\s*\(\s*per element\s*\)\s*:\s*("
                         + _NUM + r")", text):
        d["per_element"][m.group(1)] = float(m.group(2))
    d["kernel_totals"] = {}
    for m in re.finditer(r"(compute\w+)\s*\(total time\)\s*:\s*(" + _NUM
                         + r")\s*s", text):
        d["kernel_totals"][m.group(1)] = float(m.group(2))

    # IO knobs that must be OFF (these strings must be ABSENT)
    d["pgv_enabled"] = ("Peak ground motion" in text)
    d["dr_output_disabled"] = ("No dynamic rupture output enabled" in text)
    d["completed"] = (d["t_compute"] is not None)
    d["tfinal_log"] = None
    m = re.findall(r"End simulation epoch\.\s*\(at\s*(" + _NUM + r")\s*s", text)
    if m:
        d["tfinal_log"] = float(m[-1])
    return d


# ---------------------------------------------------------------------------
# Optional deck parsing (regex, not a TOML/par parser -- stdlib has no TOML
# reader before 3.11 and Expanse's system python is older).
# ---------------------------------------------------------------------------

def parse_mfem_deck(path):
    text = read_text(path)
    if text is None:
        return None
    out = {"path": path}
    for key in ("order", "ader_order", "cfl", "cfl_dg_safety",
                "checkpoint_every_steps"):
        m = re.search(r"^\s*" + key + r"\s*=\s*(" + _NUM + r")", text, re.M)
        if m:
            out[key] = float(m.group(1))
    for key in ("lts", "mixed_flux", "tag", "paraview_fault",
                "paraview_free_surface", "time_integrator"):
        m = re.search(r"^\s*" + key + r'\s*=\s*"([^"]*)"', text, re.M)
        if m:
            out[key] = m.group(1)
    return out


def parse_seissol_par(path):
    text = read_text(path)
    if text is None:
        return None
    stripped = re.sub(r"!.*", "", text)          # drop Fortran-namelist comments
    out = {"path": path}
    for key in ("ClusteredLTS", "CFL", "EndTime", "t_0", "OutputPointType",
                "SurfaceOutput", "ReceiverOutput", "EnergyOutput",
                "EnergyTerminalOutput", "Checkpoint"):
        m = re.search(r"^\s*" + key + r"\s*=\s*(" + _NUM + r")",
                      stripped, re.M | re.I)
        if m:
            out[key] = float(m.group(1))
    return out


# ---------------------------------------------------------------------------
# Leg assembly
# ---------------------------------------------------------------------------

class Leg(object):
    """One of the six measurement points."""

    def __init__(self, code, p):
        self.code = code                 # 'mfem' | 'seissol'
        self.p = p                       # MFEM polynomial degree (1/2/3)
        self.o = p + 1                   # SeisSol ConvergenceOrder (2/3/4)
        self.key = "%s_p%d" % (code, p)
        self.log = None                  # parsed log dict
        self.cali = None                 # parsed Caliper row (MFEM only)
        self.deck = None
        # (quantity, why, required) -- ALWAYS printed, never hidden.  `required`
        # False marks a quantity this code STRUCTURALLY cannot emit (it is
        # substituted from the pinned contract and gated), so it is disclosed
        # but does not by itself make the baseline incomplete.
        self.missing = []
        self.gates = []                  # (name, status, detail)
        self.m = {}                      # derived metrics

    @property
    def label(self):
        return ("MFEM p%d (ADER %d)" % (self.p, self.p + 1) if self.code == "mfem"
                else "SeisSol o%d (= p%d)" % (self.o, self.p))

    def note(self, quantity, why, required=True):
        self.missing.append((quantity, why, required))

    def gate(self, name, ok, detail=""):
        self.gates.append((name, "PASS" if ok else "FAIL", detail))
        return ok

    def warn(self, name, detail):
        self.gates.append((name, "WARN", detail))


def build_legs(run_dir, region):
    """Discover, classify and associate every artifact under run_dir."""
    legs = {}
    for code in ("mfem", "seissol"):
        for p in ORDERS:
            leg = Leg(code, p)
            legs[leg.key] = leg

    mfem_logs, seissol_logs, cali_reports = [], [], []
    decks_mfem, decks_seissol = [], []
    scanned = 0

    for path, text in scan(run_dir):
        scanned += 1
        kind = classify(path, text)
        if kind == "mfem":
            mfem_logs.append((path, parse_mfem_log(path, text)))
        elif kind == "seissol":
            seissol_logs.append((path, parse_seissol_log(path, text)))
        elif kind == "caliper":
            row = parse_caliper_report(path, text, region)
            cali_reports.append((path, row))

    for dirpath, dirnames, filenames in os.walk(run_dir):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            if fn.endswith(".toml"):
                decks_mfem.append(os.path.join(dirpath, fn))
            elif fn.endswith(".par"):
                decks_seissol.append(os.path.join(dirpath, fn))

    # --- MFEM logs: order comes from the banner, then the path token --------
    for path, d in mfem_logs:
        p = d.get("order")
        if p not in ORDERS:
            p_tok, _ = order_from_path(path)
            p = p_tok
            if p in ORDERS:
                d["notes"].append(
                    "order taken from the PATH (%s); the log has no readable "
                    "'fe order:' banner" % path)
        if p not in ORDERS:
            continue
        leg = legs["mfem_p%d" % p]
        if leg.log is None or (not leg.log.get("completed") and d.get("completed")):
            leg.log = d

    # --- SeisSol logs: exe token, then path token, then dt matching ---------
    for path, d in seissol_logs:
        o = d.get("order")
        if o not in (2, 3, 4):
            _, o_tok = order_from_path(path)
            o = o_tok
            if o in (2, 3, 4):
                d["notes"].append("order taken from the PATH token (no "
                                  "'seissol-elastic-oN-f64' string in the log)")
        if o not in (2, 3, 4) and d.get("dt"):
            for pp in ORDERS:
                if abs(d["dt"] - dt_pinned(pp)) <= 1e-3 * dt_pinned(pp):
                    o = pp + 1
                    d["notes"].append(
                        "order INFERRED from the printed dt (%.6e s) matching "
                        "the pre-registered o%d value -- no exe/path token found"
                        % (d["dt"], o))
                    break
        if o not in (2, 3, 4):
            continue
        leg = legs["seissol_p%d" % (o - 1)]
        if leg.log is None or (not leg.log.get("completed") and d.get("completed")):
            leg.log = d

    # --- Caliper reports -> MFEM legs --------------------------------------
    #  1. explicit p-token in the report path
    #  2. same stem as the leg's log (run_mfem_p1.log <-> run_mfem_p1.cali-*.txt)
    #  3. the only report in the same directory as the leg's log
    for path, row in cali_reports:
        if row is None:
            continue
        p_tok, _ = order_from_path(path)
        if p_tok in ORDERS:
            leg = legs["mfem_p%d" % p_tok]
            if leg.cali is None or "cali" in os.path.basename(path):
                leg.cali = row
            continue
        for p in ORDERS:
            leg = legs["mfem_p%d" % p]
            if leg.log is None or leg.cali is not None:
                continue
            log_path = leg.log["path"]
            stem = os.path.splitext(os.path.basename(log_path))[0]
            same_dir = os.path.dirname(path) == os.path.dirname(log_path)
            if os.path.basename(path).startswith(stem) or same_dir:
                leg.cali = row
                break

    # --- decks --------------------------------------------------------------
    for p in ORDERS:
        leg = legs["mfem_p%d" % p]
        for dk in decks_mfem:
            parsed = parse_mfem_deck(dk)
            if parsed and int(parsed.get("order", -1)) == p:
                leg.deck = parsed
                break
        leg2 = legs["seissol_p%d" % p]
        cands = decks_seissol
        if leg2.log is not None:
            same = [dk for dk in decks_seissol
                    if os.path.dirname(dk) == os.path.dirname(leg2.log["path"])]
            if same:
                cands = same
        if len(cands) == 1:
            leg2.deck = parse_seissol_par(cands[0])
        elif cands:
            tok = [dk for dk in cands if order_from_path(dk)[1] == leg2.o]
            if tok:
                leg2.deck = parse_seissol_par(tok[0])

    return legs, {"scanned": scanned, "mfem_logs": len(mfem_logs),
                  "seissol_logs": len(seissol_logs),
                  "cali_reports": len(cali_reports),
                  "mfem_decks": len(decks_mfem),
                  "seissol_decks": len(decks_seissol)}


# ---------------------------------------------------------------------------
# Metrics + gates
# ---------------------------------------------------------------------------

def conventional_paths(run_dir, leg):
    """The file(s) this leg WOULD live in under the documented convention."""
    # These are the names the two sbatches ACTUALLY write, so a "missing leg"
    # message names a path that could really have existed:
    #   run_mfem_gts_orders.sbatch      -> <OUT_ROOT>/run_p<N>.log
    #                                      <OUT_ROOT>/run_p<N>.cali-region-report.txt
    #   run_seissol_gts_orders.sbatch   -> <RUN_DIR>/run_o<N>.<jobid>.log
    # The two codes run as separate jobs on different filesystems, so gather the
    # SeisSol logs next to the MFEM ones before running this script:
    #   rsync <expanse-lustre-run-dir>/run_o*.log  <mfem runs/<jobid>>/
    if leg.code == "mfem":
        return [os.path.join(run_dir, "run_p%d.log" % leg.p),
                os.path.join(run_dir, "run_p%d.cali-region-report.txt" % leg.p)]
    return [os.path.join(run_dir, "run_o%d.<jobid>.log" % leg.o)]


def derive(leg, args):
    """Fill leg.m with the fair-metric quantities; record every gap."""
    m = leg.m
    log = leg.log
    tfinal = args.tfinal
    m["tfinal"] = tfinal
    m["dt_expected"] = dt_pinned(leg.p)
    m["steps_expected"] = steps_pinned(leg.p, tfinal)
    m["modes"] = MODES[leg.p]
    m["elements_expected"] = args.elements

    if log is None:
        leg.note("EVERYTHING", "no log for this leg")
        return

    # ---- dt --------------------------------------------------------------
    m["dt"] = log.get("dt")
    if m["dt"] is None:
        leg.note("dt", "log has no dt line "
                       "('[time] dt =' / 'Minimum timestep:')")

    # ---- steps -----------------------------------------------------------
    if m["dt"]:
        m["steps_derived"] = int(math.ceil(tfinal / m["dt"]))
    else:
        m["steps_derived"] = None

    if leg.code == "mfem":
        m["steps_progress"] = log.get("steps_progress")
        m["steps_planned"] = log.get("nsteps_planned")
        m["steps_caliper"] = None
        if leg.cali and leg.cali.get("calls") is not None and log.get("ranks"):
            raw = leg.cali["calls"]
            per_rank = raw / float(log["ranks"])
            tgt = m["steps_derived"] or m["steps_expected"]
            if abs(raw - tgt) <= max(2, 0.01 * tgt):
                m["steps_caliper"] = raw
                m["caliper_calls_scope"] = "already per-rank"
            elif abs(per_rank - tgt) <= max(2, 0.01 * tgt):
                m["steps_caliper"] = int(round(per_rank))
                m["caliper_calls_scope"] = ("rank-summed (raw %d / %d ranks)"
                                            % (raw, log["ranks"]))
            else:
                m["caliper_calls_scope"] = (
                    "UNRECOGNISED: raw=%d, /ranks=%.1f, expected~%s"
                    % (raw, per_rank, tgt))
                leg.note("steps (Caliper region.count)",
                         "Calls=%d matches neither the derived step count nor "
                         "steps x ranks; the step count fell back to the "
                         "'step k/N' progress line" % raw, required=False)
        m["steps"] = (m["steps_caliper"] or m["steps_progress"]
                      or m["steps_derived"])
        if m["steps"] is None:
            leg.note("steps", "no Caliper Calls, no 'step k/N' progress line, "
                              "and no dt to derive from")
    else:
        m["steps"] = m["steps_derived"]
        lp = log.get("last_progress")
        m["steps_last_progress"] = lp[0] if lp else None
        if m["steps"] is None:
            leg.note("steps", "no 'Minimum timestep:' line -> cannot derive "
                              "ceil(tfinal/dt); SeisSol never prints a total")

    # ---- elements --------------------------------------------------------
    if leg.code == "seissol":
        m["elements"] = log.get("elements") or args.elements
        if log.get("elements") is None:
            leg.note("elements (measured)",
                     "no cluster histogram / 'Cell count:' line; using the "
                     "pinned %d" % args.elements, required=False)
    else:
        # The MFEM driver never prints a GLOBAL element count (it prints local
        # NE only under --lts-report, spatial_dyn_driver.cpp:3504).  The mesh is
        # md5-identical to SeisSol's, so the pinned count is used and gated
        # against SeisSol's parsed value below.
        m["elements"] = args.elements
        leg.note("elements (measured)",
                 "seas_spatial_dyn_driver does not print a global element "
                 "count (local NE only, and only under --lts-report, "
                 "spatial_dyn_driver.cpp:3504); using the pinned/md5-verified "
                 "%d, cross-gated against SeisSol's parsed 'Cell count:'"
                 % args.elements, required=False)

    # ---- cores -----------------------------------------------------------
    if leg.code == "mfem":
        ranks = log.get("ranks")
        if args.mfem_cores:
            m["cores_alloc"] = m["cores_active"] = args.mfem_cores
            m["cores_src"] = "--mfem-cores override"
        elif ranks:
            # pure MPI, --cpus-per-task=1 -> 1 core per rank, none idle
            m["cores_alloc"] = m["cores_active"] = ranks
            m["cores_src"] = "'ranks:' banner (pure MPI, 1 core/rank)"
        else:
            m["cores_alloc"] = m["cores_active"] = None
            m["cores_src"] = None
            leg.note("cores", "no 'ranks:' banner and no --mfem-cores")
    else:
        ranks = log.get("ranks")
        thr = log.get("omp_threads")
        cpr = log.get("slurm_cpus_per_rank") or args.seissol_cpus_per_rank
        if ranks and thr:
            m["cores_active"] = ranks * thr
            m["cores_alloc"] = (args.seissol_cores if args.seissol_cores
                                else ranks * cpr)
            m["cores_src"] = ("%d ranks x %d OMP threads (active); x %d "
                              "cpus/rank (allocated)" % (ranks, thr, cpr))
        else:
            m["cores_active"] = None
            m["cores_alloc"] = args.seissol_cores
            m["cores_src"] = "--seissol-cores override" if args.seissol_cores else None
            leg.note("cores", "no 'Using MPI with #ranks:' / 'Using OpenMP "
                              "with #threads/rank:' lines")

    # ---- T_loop ----------------------------------------------------------
    if leg.code == "mfem":
        if leg.cali is None:
            leg.note("T_loop", "no Caliper region report for '%s'.  "
                               "MFEM_PERF_SCOPE compiles to NOTHING without "
                               "libcaliper (general/annotation.hpp:31) and the "
                               "driver has no other timer -- rebuild with "
                               "USE_CALIPER=YES bash build_expanse.sh"
                     % DEFAULT_REGION)
            m["t_loop"] = None
        elif leg.cali.get("max") is None:
            leg.note("T_loop", "Caliper row found but no time column could be "
                               "mapped (row: %s)" % leg.cali.get("row"))
            m["t_loop"] = None
        else:
            m["t_loop"] = leg.cali["max"]
            m["t_loop_mean"] = leg.cali.get("avg")
            m["t_loop_stat"] = "max"
            m["t_loop_src"] = ("Caliper '%s' %s Max-across-ranks (%s mapping)"
                               % (DEFAULT_REGION,
                                  "inclusive" if leg.cali["inclusive"] else "exclusive",
                                  leg.cali["mapping_mode"]))
            if leg.cali.get("avg"):
                m["rank_imbalance"] = leg.cali["max"] / leg.cali["avg"]
    else:
        tc = log.get("t_compute")
        if tc is None:
            leg.note("T_loop", "no terminal 'Simulation time (compute):' line "
                               "(Simulator.cpp:150) -- run did not finish "
                               "cleanly?  The per-epoch '(compute)' line is "
                               "NOT a substitute: it double-counts IO.")
            m["t_loop"] = None
        else:
            # HEADLINE on the MAX, to match MFEM's Caliper Max-across-ranks.
            # Stopwatch::print LEADS WITH THE MEAN (Monitoring/Stopwatch.cpp:69-76),
            # and quoting SeisSol-mean against MFEM-max would credit SeisSol with
            # its whole load imbalance (~6.4% here) -- the same magnitude as the
            # 240-vs-256 core caveat.  Both are kept so section 4 can print the
            # ratio max-vs-max AND mean-vs-mean.
            if tc.get("max"):
                m["t_loop"] = tc["max"]
                m["t_loop_stat"] = "max"
                m["t_loop_src"] = ("terminal 'Simulation time (compute):' per-rank "
                                   "MAX (from the 'range: [min, max]' field), chosen "
                                   "to match MFEM's Caliper Max-across-ranks")
            else:
                m["t_loop"] = tc["mean"]
                m["t_loop_stat"] = "mean"
                m["t_loop_src"] = ("terminal 'Simulation time (compute):' per-rank "
                                   "MEAN -- no 'range: [min, max]' field was parsed, "
                                   "so this is NOT statistic-matched to MFEM's Max")
            m["t_loop_mean"] = tc.get("mean")
            if tc.get("max") and tc.get("mean"):
                m["rank_imbalance"] = tc["max"] / tc["mean"]

    # ---- the metrics -----------------------------------------------------
    if m.get("t_loop") and tfinal > 0:
        m["s_per_sim_s"] = m["t_loop"] / tfinal
        m["sim_s_per_wall_h"] = 3600.0 / m["s_per_sim_s"]
    else:
        m["s_per_sim_s"] = None
        m["sim_s_per_wall_h"] = None
        leg.note("s/sim-s and sim-s/wall-h", "T_loop unavailable")

    if m.get("steps") and m.get("elements"):
        m["updates"] = m["steps"] * m["elements"]
    else:
        m["updates"] = None
        leg.note("element-updates U = S x E", "steps and/or elements missing")

    for tag, ckey in (("us_core_alloc", "cores_alloc"),
                      ("us_core_active", "cores_active")):
        if m.get("t_loop") and m.get("updates") and m.get(ckey):
            m[tag] = m["t_loop"] * m[ckey] * 1e6 / m["updates"]
            m[tag.replace("us_core", "ns_core_dof")] = (
                m[tag] * 1000.0 / m["modes"])
        else:
            m[tag] = None
            if m.get(ckey) is None and m.get("t_loop") and m.get("updates"):
                leg.note("us-core/update (%s)" % ckey, "core count unknown")

    if m.get("t_loop") and m.get("steps"):
        m["s_per_step"] = m["t_loop"] / m["steps"]


def check_gates(leg, args, seissol_elements):
    """Fairness gates.  Every failure is printed with its evidence."""
    log = leg.log
    if log is None:
        leg.gate("leg present", False, "no log found")
        return

    tol = args.dt_tol
    dt_exp = leg.m["dt_expected"]
    dt = leg.m.get("dt")

    # G1 -- dt matches the pre-registered CFL prediction  (R1 / R3)
    if dt is None:
        leg.gate("G1 dt == pre-registered %.6e s" % dt_exp, False,
                 "dt not parsed")
    else:
        ok = abs(dt - dt_exp) <= tol * dt_exp
        leg.gate("G1 dt == pre-registered %.6e s" % dt_exp, ok,
                 "measured %.6e s (rel dev %.2e)" % (dt, abs(dt - dt_exp) / dt_exp))

    # G2 -- steps == ceil(tfinal/dt)
    s, sd = leg.m.get("steps"), leg.m.get("steps_derived")
    if s and sd:
        ok = abs(s - sd) <= max(1, 0.005 * sd)
        leg.gate("G2 steps == ceil(tfinal/dt)", ok,
                 "measured %d vs derived %d" % (s, sd))
    else:
        leg.gate("G2 steps == ceil(tfinal/dt)", False, "steps unavailable")

    # G3 -- elements identical to the other code / the pin  (R1)
    e = leg.m.get("elements")
    ok = (e == args.elements)
    detail = "%s vs pinned %d" % (e, args.elements)
    if seissol_elements and e != seissol_elements:
        ok = False
        detail += "; SeisSol parsed %d" % seissol_elements
    leg.gate("G3 elements == %d (same mesh)" % args.elements, ok, detail)

    if leg.code == "mfem":
        # G4 -- accuracy-matched pairing p <-> ader p+1  (R2)
        ao = log.get("ader_order")
        leg.gate("G4 ader_order == p+1 = %d" % (leg.p + 1),
                 ao == leg.p + 1, "banner 'ader order: %s'" % ao)
        # G5 -- cfl_dg_safety == 1.0 (else dt/3 and a 3x-inflated MFEM number)
        sf = log.get("cfl_dg_safety")
        leg.gate("G5 cfl_dg_safety == 1.0", sf is not None and abs(sf - 1.0) < 1e-9,
                 "log reports cfl_dg_safety = %s (deck default is 3.0!)" % sf)
        # G6 -- Courant number
        leg.gate("G6 cfl == 0.5", log.get("cfl") is not None
                 and abs(log["cfl"] - 0.5) < 1e-9, "banner 'cfl: %s'" % log.get("cfl"))
        # G7 -- GTS
        leg.gate("G7 GTS (no '[lts] ENABLED' banner)",
                 not log.get("lts_enabled_banner"),
                 "'[lts] ENABLED' %s" % ("PRESENT" if log.get("lts_enabled_banner")
                                         else "absent"))
        # G8 -- IO off inside the timed region  (R4).  TWO conditions: the
        # ParaView master gate AND the free-surface slice, which sits outside it.
        leg.gate("G8a ParaView OFF", bool(log.get("paraview_off")),
                 "'ParaView output: OFF' %s" %
                 ("found" if log.get("paraview_off") else "NOT found -- writes "
                  "are inside the timed step region (spatial_dyn_driver.cpp:5121)"))
        fs_deck = leg.deck.get("paraview_free_surface") if leg.deck else None
        fs_ok = (not log.get("free_surface_banner")) and (fs_deck in (None, "off"))
        leg.gate("G8b free-surface slice OFF", fs_ok,
                 "log '  FreeSurf:' banner %s; deck paraview_free_surface = %s "
                 "(NOT under the paraview_enabled master gate; DEFAULT is \"vtu\" "
                 "@0.05 s, written from inside seas::spatial_dyn::step)"
                 % ("PRESENT" if log.get("free_surface_banner") else "absent",
                    fs_deck if fs_deck is not None else "<deck not parsed>"))
        # G9 -- no IN-LOOP checkpoint.  The old log-string test was UNREACHABLE on
        # the GTS path: '[checkpoint] ... written' is printed only by the LTS
        # branch (spatial_dyn_driver.cpp:5258); the GTS in-loop write at
        # :5209-5214 is completely silent, so that gate always PASSed even for a
        # leg that fired a full-state 256-rank write inside the timed region.
        # Gate on the DECK value instead, which is the thing that decides it.
        ces = leg.deck.get("checkpoint_every_steps") if leg.deck else None
        steps = leg.m.get("steps")
        if ces is None:
            leg.warn("G9 no in-loop checkpoint",
                     "deck checkpoint_every_steps not parsed -- cannot verify that "
                     "no in-loop checkpoint fired (the log gives no signal on the "
                     "GTS path).  Count cp_checkpoint_r*.txt: exactly <nranks> is "
                     "the expected post-loop final write; more means an in-loop one.")
        else:
            ok = (steps is None) or (ces > steps)
            leg.gate("G9 no in-loop checkpoint", ok,
                     "deck checkpoint_every_steps = %g vs %s steps -- %s"
                     % (ces, steps if steps is not None else "?",
                        "no in-loop write can fire" if ok else
                        "an in-loop full-state write DID fire inside the timed region"))
        if log.get("checkpoint_written"):
            leg.warn("checkpoint (LTS signal)",
                     "'[checkpoint] ... written' present -- that string is emitted "
                     "only by the LTS branch, so this leg may not be GTS.")
        # R4 SYMMETRY: a non-empty [problem].tag wires a per-step station writer
        # INSIDE the timed scope while the SeisSol legs write nothing at all.
        tagv = leg.deck.get("tag") if leg.deck else None
        if tagv:
            leg.gate("G9b [problem].tag empty (IO symmetry)", False,
                     "deck sets tag = '%s' -> TPV station writer runs every step "
                     "inside seas::spatial_dyn::step (spatial_dyn_driver.cpp:4467 "
                     "-> :5126), fflushing 9 files on <=9 ranks, with a row count "
                     "that GROWS with order.  The SeisSol partner legs run "
                     "OutputPointType = 0 / ReceiverOutput = 0 -- zero writes.  "
                     "This leg is NOT IO-comparable to its SeisSol row." % tagv)
        elif leg.deck:
            leg.gate("G9b [problem].tag empty (IO symmetry)", True,
                     "tag is empty -> no station writer wired")
        # time_integrator: the RK path uses RkCflFactor = 3/(2p+1), i.e. dt x3.
        ti = leg.deck.get("time_integrator") if leg.deck else None
        if ti is not None and ti != "ader":
            leg.gate("G9c time_integrator == 'ader'", False,
                     "deck sets time_integrator = '%s'; the RK branch takes dt x3 "
                     "(spatial_friction.hpp:767 vs :744) -- 3x FEWER steps than "
                     "SeisSol, while 'cfl_dg_safety = 1' still prints." % ti)
        if log.get("mixed_flux") not in (None, "none"):
            leg.warn("mixed_flux", "expected 'none' under ADER, log says '%s' "
                                   "(ComputeMaxDt ABORTS on 'adjacent' under "
                                   "ADER, wave_operator.inl:7717-7721)"
                     % log.get("mixed_flux"))
        if leg.deck:
            if str(leg.deck.get("lts", "off")) != "off":
                leg.gate("G7b deck [numerics].lts == 'off'", False,
                         "deck says lts = '%s'" % leg.deck.get("lts"))
            if leg.deck.get("tag"):
                leg.warn("per-step SCEC station writes",
                         "deck [problem].tag = '%s' wires the station writer, "
                         "which fflushes every step on <=9 ranks INSIDE the "
                         "timed region (spatial_dyn_driver.cpp:4468, :5126); "
                         "declare this asymmetry (p3 writes %d rows/station, "
                         "p1 only %d)" % (leg.deck["tag"],
                                          steps_pinned(3, args.tfinal),
                                          steps_pinned(1, args.tfinal)))
        if log.get("V_max_peak"):
            leg.warn("physics comparability",
                     "V_max peak %.4g m/s (p3 reference ~17.3 m/s at t~1.48 s). "
                     "A leg that does not rupture is cheap for a PHYSICS reason "
                     "and its throughput is not comparable."
                     % log["V_max_peak"])
        else:
            leg.warn("physics comparability",
                     "no 'V_max' progress line parsed -- cannot confirm this leg "
                     "ran the same rupture as the others")
    else:
        # G4 -- GTS, three independent ways  (R1)
        leg.gate("G4 'GTS has been selected.'", bool(log.get("gts_selected")),
                 "LtsWeights.cpp:124")
        cl = log.get("speedup_clts")
        leg.gate("G5 clustered-LTS speedup == 1", cl is not None and abs(cl - 1.0) < 1e-6,
                 "'Theoretical speedup to GTS: ... %s clustered LTS'" % cl)
        nc = log.get("n_clusters")
        leg.gate("G6 single cluster-histogram row", nc == 1,
                 "%s row(s) -- >1 means LTS is still on" % nc)
        # G7 -- exactly one epoch  (R4: output cadence changes the STEP COUNT)
        ne = log.get("epochs_start")
        leg.gate("G7 exactly ONE simulation epoch", ne == 1,
                 "%s 'Start simulation epoch.' line(s); >1 means a periodic "
                 "output module survived and clipped dt at every sync point"
                 % ne)
        # G8 -- IO off
        leg.gate("G8 PGV sampling OFF", not log.get("pgv_enabled"),
                 "'Peak ground motion (PGV) enabled' %s -- it folds every free-"
                 "surface face at EVERY local timestep INSIDE the compute "
                 "stopwatch (TimeCluster.cpp:905)"
                 % ("PRESENT" if log.get("pgv_enabled") else "absent"))
        io = log.get("t_io")
        if io and io.get("mean") is not None:
            frac = None
            tot = log.get("t_total")
            if tot and tot.get("mean"):
                frac = 100.0 * io["mean"] / tot["mean"]
            leg.gate("G9 blocking IO negligible",
                     frac is not None and frac < 1.0,
                     "blocking IO %.3f s%s" % (io["mean"],
                                               "" if frac is None else
                                               " = %.2f %% of total" % frac))
        # audit trail on the derived step count
        lp = leg.m.get("steps_last_progress")
        if lp is not None and leg.m.get("steps"):
            leg.gate("G10 last progress line within 100 of derived steps",
                     abs(lp - leg.m["steps"]) < 100,
                     "last 'Max cluster / LTS cycle updates since sync: %d' vs "
                     "derived %d" % (lp, leg.m["steps"]))
        else:
            leg.warn("G10 step audit trail",
                     "no 'Max cluster / LTS cycle updates since sync:' line "
                     "parsed -- the derived step count is unaudited")
        if leg.deck:
            if leg.deck.get("ClusteredLTS") is not None:
                leg.gate("G4b deck ClusteredLTS == 1",
                         abs(leg.deck["ClusteredLTS"] - 1.0) < 1e-9,
                         "par says %s" % leg.deck["ClusteredLTS"])
            if leg.deck.get("CFL") is not None:
                leg.gate("G11 deck CFL == 0.5",
                         abs(leg.deck["CFL"] - 0.5) < 1e-9,
                         "par says %s" % leg.deck["CFL"])
            for key in ("OutputPointType", "SurfaceOutput", "ReceiverOutput",
                        "EnergyOutput", "EnergyTerminalOutput", "Checkpoint"):
                if leg.deck.get(key) is not None and leg.deck[key] != 0.0:
                    leg.gate("G12 %s == 0" % key, False,
                             "par says %s -- this is INSIDE the compute "
                             "stopwatch (EnergyOutput additionally changes the "
                             "DR friction KERNEL, DRParameters.cpp:154)"
                             % leg.deck[key])
        if log.get("tfinal_log") is not None:
            leg.gate("G13 reached EndTime %.3f s" % args.tfinal,
                     abs(log["tfinal_log"] - args.tfinal) < 1e-6,
                     "last epoch ended at %s s" % log["tfinal_log"])


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def fmt(v, spec="%.4g", na="  --  "):
    if v is None:
        return na
    try:
        return spec % v
    except (TypeError, ValueError):
        return str(v)


def rule(ch="-", n=100):
    return ch * n


def print_contract(args):
    print(rule("="))
    print("0. PRE-REGISTERED CONTRACT  (fixed BEFORE any timing number is read)")
    print(rule("="))
    print("  Mesh          : SCEC TPV104 200 m, %d tets (md5-identical .msh / .puml.h5)"
          % args.elements)
    print("  Window        : t = 0 .. %.3f s        Courant : 0.5 (both codes)" % args.tfinal)
    print("  Scheme        : GTS both sides  (MFEM [numerics].lts=\"off\"; SeisSol ClusteredLTS=1)")
    print("  Pairing (R2)  : ACCURACY-MATCHED  p_k  <->  o_(k+1)   (SeisSol o = p + 1)")
    print("                  At the matched order the two codes advance state vectors of")
    print("                  IDENTICAL size (%d/%d/%d modes x 9 quantities per element) with"
          % (MODES[1], MODES[2], MODES[3]))
    print("                  IDENTICAL dt and step counts on IDENTICAL geometry.")
    print("  dt (R3)       : NOT held equal across orders -- each order gets its own")
    print("                  CFL-derived dt; both codes use the same algebraic formula,")
    print("                  dt propto 1/(2p+1) == 1/(2o-1).")
    print("")
    hdr = ("  %-6s %-5s %-6s %-6s %-6s %14s %8s %16s" %
           ("MFEM", "SeisSol", "modes", "faultQP", "sub", "dt [s]", "steps",
            "element-updates"))
    print(hdr)
    print("  " + rule("-", len(hdr) - 2))
    for p in ORDERS:
        print("  p%-5d o%-4d %-6d %-7d %-6d %14.6e %8d %16.5e"
              % (p, p + 1, MODES[p], FAULT_QP_PER_FACE[p], p + 1,
                 dt_pinned(p), steps_pinned(p, args.tfinal),
                 float(steps_pinned(p, args.tfinal) * args.elements)))
    print("")
    print("  'sub' = friction sub-steps per macro step (MFEM: O = ader_order,")
    print("          spatial_dyn_driver.cpp:4189-4192; SeisSol: ConvergenceOrder DR")
    print("          time-quadrature points).  Identical under the matched pairing.")
    print("  dt anchor: the verified p3/o4 dt = %.6e s at cfl 0.5; p1,p2 = x7/3, x7/5."
          % DT_REF_P3)


def print_inventory(run_dir, legs, stats):
    print("")
    print(rule("="))
    print("1. LEG INVENTORY")
    print(rule("="))
    print("  run dir: %s" % os.path.abspath(run_dir))
    print("  scanned %d text artifact(s): %d MFEM log(s), %d SeisSol log(s), "
          "%d Caliper report(s), %d TOML deck(s), %d .par deck(s)"
          % (stats["scanned"], stats["mfem_logs"], stats["seissol_logs"],
             stats["cali_reports"], stats["mfem_decks"], stats["seissol_decks"]))
    print("")
    missing_any = False
    for code in ("mfem", "seissol"):
        for p in ORDERS:
            leg = legs["%s_p%d" % (code, p)]
            if leg.log is None:
                missing_any = True
                print("  [MISSING] %-22s  NO LOG FOUND." % leg.label)
                for cp in conventional_paths(run_dir, leg):
                    print("            expected (convention): %s" % cp)
                print("            (discovery is content-based: an MFEM log must contain")
                print("             'fe order:' or '[time] dt ='; a SeisSol log must contain")
                print("             'Minimum timestep:' or 'Cell count:'.)")
                continue
            extra = ""
            if leg.code == "mfem":
                if leg.cali is None:
                    missing_any = True
                    extra = "   [MISSING Caliper region report -> NO T_loop]"
                else:
                    extra = "   + %s" % os.path.relpath(leg.cali["path"], run_dir)
            print("  [found]   %-22s  %s%s"
                  % (leg.label, os.path.relpath(leg.log["path"], run_dir), extra))
            if leg.deck:
                print("            deck: %s" % os.path.relpath(leg.deck["path"], run_dir))
            for n in leg.log.get("notes", []):
                print("            note: %s" % n)
            if leg.code == "mfem" and leg.cali is None:
                for cp in conventional_paths(run_dir, leg)[1:]:
                    print("            expected (convention): %s" % cp)
    return missing_any


def print_gates(legs):
    print("")
    print(rule("="))
    print("2. FAIRNESS GATES  (R1-R6)   -- a FAIL voids that leg's place in the baseline")
    print(rule("="))
    n_fail = 0
    for code in ("mfem", "seissol"):
        for p in ORDERS:
            leg = legs["%s_p%d" % (code, p)]
            print("")
            print("  %s" % leg.label)
            if not leg.gates:
                print("    (no gates evaluated -- leg absent)")
                continue
            for name, status, detail in leg.gates:
                if status == "FAIL":
                    n_fail += 1
                print("    [%-4s] %-46s %s" % (status, name, detail))
    return n_fail


def print_baseline(legs, args):
    print("")
    print(rule("="))
    print("3. BASELINE TABLE   (T_loop = time-integration region ONLY; setup excluded)")
    print(rule("="))
    hdr = ("  %-22s %12s %7s %11s %11s %12s %11s %11s %11s"
           % ("leg", "dt [s]", "steps", "T_loop [s]", "s / sim-s",
              "sim-s/wall-h", "cores(a/A)", "us-core/upd", "us-core/upd"))
    print(hdr)
    print("  %-22s %12s %7s %11s %11s %12s %11s %11s %11s"
          % ("", "", "", "", "", "", "alloc/active", "@C_alloc", "@C_active"))
    print("  " + rule("-", len(hdr) - 2))
    for code in ("mfem", "seissol"):
        for p in ORDERS:
            leg = legs["%s_p%d" % (code, p)]
            m = leg.m
            cores = "%s/%s" % (fmt(m.get("cores_alloc"), "%d", "?"),
                               fmt(m.get("cores_active"), "%d", "?"))
            print("  %-22s %12s %7s %11s %11s %12s %11s %11s %11s"
                  % (leg.label,
                     fmt(m.get("dt"), "%.6e"),
                     fmt(m.get("steps"), "%d"),
                     fmt(m.get("t_loop"), "%.1f"),
                     fmt(m.get("s_per_sim_s"), "%.1f"),
                     fmt(m.get("sim_s_per_wall_h"), "%.4f"),
                     cores,
                     fmt(m.get("us_core_alloc"), "%.2f"),
                     fmt(m.get("us_core_active"), "%.2f")))
        if code == "mfem":
            print("  " + rule("-", len(hdr) - 2))
    print("")
    print("  us-core/upd = T_loop * cores * 1e6 / (steps * elements)   [R6]")
    print("  C_alloc  = cores the job BILLED (billing-fair headline, identical shape-wise)")
    print("  C_active = cores actually executing.  The validated SeisSol shape sets")
    print("             OMP_NUM_THREADS = cpus_per_task - 1 AND SEISSOL_COMMTHREAD=0, so")
    print("             one core per rank is allocated but idle (16x15 = 240 of 256).")
    print("             Do NOT 'fix' this by raising OMP_NUM_THREADS -- R5 forbids changing")
    print("             either validated shape; the idle core is itself a documented finding.")
    print("")
    print("")
    print("  PROVENANCE of every T_loop / steps / cores cell above (R6 audit trail):")
    for code in ("mfem", "seissol"):
        for p in ORDERS:
            leg = legs["%s_p%d" % (code, p)]
            m = leg.m
            if leg.log is None:
                print("    %-22s leg absent" % leg.label)
                continue
            print("    %-22s T_loop: %s%s"
                  % (leg.label, m.get("t_loop_src", "UNAVAILABLE"),
                     ("   [rank imbalance Max/Avg = %.3f]" % m["rank_imbalance"])
                     if m.get("rank_imbalance") else ""))
            if leg.code == "mfem":
                print("    %-22s   steps: caliper=%s, progress-line=%s, planned=%s%s"
                      % ("", m.get("steps_caliper"), m.get("steps_progress"),
                         m.get("steps_planned"),
                         ("; region.count is %s" % m["caliper_calls_scope"])
                         if m.get("caliper_calls_scope") else ""))
            else:
                print("    %-22s   steps: derived ceil(tfinal/dt)=%s, last progress line=%s"
                      % ("", m.get("steps_derived"), m.get("steps_last_progress")))
            print("    %-22s   cores: %s" % ("", m.get("cores_src") or "UNKNOWN"))
    print("")
    print("    rank imbalance = Max/Avg across ranks.  Well above 1 means the quoted wall is")
    print("    one straggler rank rather than the whole job -- disclose it, do not average it")
    print("    away (MFEM uses Max because a bulk-synchronous step ends when the last rank")
    print("    arrives; SeisSol's own figure is a per-rank MEAN, so the two are not identically")
    print("    conservative and the imbalance column is how the reader sees the difference).")
    print("")
    print("  per-DOF-update view (us-core/upd divided by modes/element -> ns-core per")
    print("  mode-update; within-code lens only -- at a matched order the two codes have")
    print("  the SAME modes/element, so this adds no cross-code information):")
    print("    %-22s %14s %16s %16s" % ("leg", "modes/elem", "ns-core @C_alloc",
                                        "ns-core @C_active"))
    for code in ("mfem", "seissol"):
        for p in ORDERS:
            leg = legs["%s_p%d" % (code, p)]
            print("    %-22s %14d %16s %16s"
                  % (leg.label, leg.m.get("modes", 0),
                     fmt(leg.m.get("ns_core_dof_alloc"), "%.2f"),
                     fmt(leg.m.get("ns_core_dof_active"), "%.2f")))


def print_cross_code(legs, args):
    print("")
    print(rule("="))
    print("4. CROSS-CODE RATIO AT MATCHED ORDER   MFEM p_k  vs  SeisSol o_(k+1)")
    print(rule("="))
    print("  STATISTIC MATCHING (this is a fairness requirement, not a detail):")
    print("    MFEM's wall is a Caliper MAX over 256 ranks; SeisSol's 'Simulation time")
    print("    (compute):' LEADS WITH THE MEAN over 16 ranks and carries the max only in its")
    print("    'range: [min, max]' field.  Quoting MFEM-max against SeisSol-mean would credit")
    print("    SeisSol with its whole load imbalance (~6.4% on this problem) -- the same")
    print("    magnitude as the 240-vs-256 core caveat.  Both pairings are printed below.")
    print("")
    hdr = ("  %-14s %14s %14s %12s %14s %14s %12s"
           % ("matched order", "MFEM s/sim-s", "SS s/sim-s", "ratio",
              "MFEM us/upd", "SS us/upd", "ratio"))
    print("  --- HEADLINE: max-vs-max (a bulk-synchronous step ends when the LAST rank arrives)")
    print(hdr)
    print("  " + rule("-", len(hdr) - 2))
    for p in ORDERS:
        a = legs["mfem_p%d" % p].m
        b = legs["seissol_p%d" % p].m
        for k, lg in ((a, legs["mfem_p%d" % p]), (b, legs["seissol_p%d" % p])):
            if k.get("t_loop") is not None and k.get("t_loop_stat") != "max":
                print("      ^ NOTE: %s contributes a %s, not a max -- this row is NOT"
                      % (lg.label, k.get("t_loop_stat")))
                print("        statistic-matched; see its provenance block in section 3.")
        r_t = (a["s_per_sim_s"] / b["s_per_sim_s"]
               if a.get("s_per_sim_s") and b.get("s_per_sim_s") else None)
        r_u = (a["us_core_alloc"] / b["us_core_alloc"]
               if a.get("us_core_alloc") and b.get("us_core_alloc") else None)
        print("  %-14s %14s %14s %12s %14s %14s %12s"
              % ("p%d <-> o%d" % (p, p + 1),
                 fmt(a.get("s_per_sim_s"), "%.1f"),
                 fmt(b.get("s_per_sim_s"), "%.1f"),
                 fmt(r_t, "%.2fx"),
                 fmt(a.get("us_core_alloc"), "%.2f"),
                 fmt(b.get("us_core_alloc"), "%.2f"),
                 fmt(r_u, "%.2fx")))
        # consistency: with equal dt, equal steps, equal elements and equal
        # C_alloc, the two ratios MUST coincide.  A gap means one of those
        # equalities silently broke.
        if r_t and r_u and abs(r_t - r_u) > 0.01 * max(r_t, r_u):
            print("      ^ WARNING: the wall-ratio and the per-update ratio disagree by "
                  "%.1f %%." % (100.0 * abs(r_t - r_u) / max(r_t, r_u)))
            print("        With equal dt, equal steps, equal elements and equal C_alloc "
                  "they must be identical;")
            print("        one of those equalities is broken -- inspect the gates above "
                  "before quoting either number.")
    print("")
    print("  --- COMPANION: mean-vs-mean (both codes' per-rank average)")
    print("  %-14s %14s %14s %12s" % ("matched order", "MFEM T_mean", "SS T_mean", "ratio"))
    print("  " + rule("-", 56))
    tf = args.tfinal if args.tfinal > 0 else 1.0
    for p in ORDERS:
        a = legs["mfem_p%d" % p].m
        b = legs["seissol_p%d" % p].m
        am, bm = a.get("t_loop_mean"), b.get("t_loop_mean")
        r = (am / bm) if (am and bm) else None
        print("  %-14s %14s %14s %12s"
              % ("p%d <-> o%d" % (p, p + 1),
                 fmt(am, "%.1f"), fmt(bm, "%.1f"), fmt(r, "%.2fx")))
    print("    (T in seconds over the %.3g s window; divide by the window for s/sim-s.)" % tf)
    print("    If the max-vs-max and mean-vs-mean ratios differ materially, the gap IS the")
    print("    two codes' load-imbalance difference -- report it, do not pick the flattering one.")
    print("")
    print("  Ratios use C_alloc (256 cores both codes).  Known MFEM-side confounds, to be")
    print("  restated wherever these ratios are quoted:")
    print("   * pmesh.SetCurvature(order) is applied UNCONDITIONALLY to a straight-sided")
    print("     tet mesh (drivers/spatial_dyn_driver.cpp:1733), so MFEM evaluates a 4/10/20-")
    print("     node H1 nodal transformation per QP where SeisSol is always affine.  Part of")
    print("     the p3 gap is a redundant GEOMETRY order, not DG arithmetic.  --face-cache")
    print("     removes it on non-fault interior faces only (wave_operator.inl:884-905).")
    print("   * instrumentation asymmetry: MFEM's timed region contains ~17 nested Caliper")
    print("     regions; SeisSol's compute stopwatch (Simulator.cpp:106-108) is a bare")
    print("     start/pause pair with no annotation at all.  This job's CALI_CONFIG omits")
    print("     profile.mpi for that reason, but the residual annotation tax is one-sided.")
    print("   * MFEM-only per-step global syncs with no SeisSol counterpart: a Q.Norml2() +")
    print("     1-value MPI_Allreduce NaN tripwire every step across all 256 ranks")
    print("     (spatial_dyn_driver.cpp:5190-5199, NOT disableable).  The order-dependent")
    print("     R-101 shared-fault check (:5155-5162, 12/20/28 firings at p1/p2/p3) is")
    print("     removed uniformly by SEAS_R101_SKIP=1 in the run sbatch.")
    print("   * partition-count asymmetry: MFEM is decomposed 256 ways, SeisSol 16 ways")
    print("     (~2.5x more halo for MFEM, all of it MPI).  Per-THREAD work is matched to")
    print("     within 6.7 % (9,628 vs 10,270 elements); only the halo differs.  The")
    print("     penalty is largest at p1 -- do not lead with the p1 ratio.")
    print("   * SeisSol's o2 build pads 4 modes to the machine vector width, so its p1 leg")
    print("     does far more HW-FLOP than NZ-FLOP; check the HW/NZ ratio below.")


def print_within_code(legs, args):
    print("")
    print(rule("="))
    print("5. WITHIN-CODE ORDER SCALING  (p1 -> p2 -> p3, each code against its own p1)")
    print(rule("="))
    print("  The total-cost ratio contains an EXACTLY KNOWN step-count penalty:")
    print("  S propto (2p+1) = 3 : 5 : 7.  It is divided out below so per-step work is")
    print("  legible on its own -- only the per-step column is informative.")
    print("")
    for code in ("mfem", "seissol"):
        base = legs["%s_p1" % code].m
        print("  %s" % ("MFEM" if code == "mfem" else "SeisSol"))
        hdr = ("    %-10s %13s %10s %13s %10s %14s %10s"
               % ("order", "s/sim-s", "vs p1", "s/step", "vs p1",
                  "us-core/upd", "vs p1"))
        print(hdr)
        print("    " + rule("-", len(hdr) - 4))
        for p in ORDERS:
            m = legs["%s_p%d" % (code, p)].m
            r_tot = (m["s_per_sim_s"] / base["s_per_sim_s"]
                     if m.get("s_per_sim_s") and base.get("s_per_sim_s") else None)
            r_stp = (m["s_per_step"] / base["s_per_step"]
                     if m.get("s_per_step") and base.get("s_per_step") else None)
            r_upd = (m["us_core_alloc"] / base["us_core_alloc"]
                     if m.get("us_core_alloc") and base.get("us_core_alloc") else None)
            print("    %-10s %13s %10s %13s %10s %14s %10s"
                  % ("p%d" % p if code == "mfem" else "o%d" % (p + 1),
                     fmt(m.get("s_per_sim_s"), "%.1f"), fmt(r_tot, "%.2fx"),
                     fmt(m.get("s_per_step"), "%.4f"), fmt(r_stp, "%.2fx"),
                     fmt(m.get("us_core_alloc"), "%.2f"), fmt(r_upd, "%.2fx")))
        print("")
    ref = [(2 * p + 1) * MODES[p] for p in ORDERS]
    print("  Reference MODEL (NOT a prediction): a cost linear in DOFs would give total")
    print("  (2p+1)*modes = %d : %d : %d  =  1.00x : %.2fx : %.2fx, and per-step modes"
          % (ref[0], ref[1], ref[2], ref[1] / float(ref[0]), ref[2] / float(ref[0])))
    print("  = 4 : 10 : 20 = 1.00x : 2.50x : 5.00x.  ADER-DG volume kernels are closer to")
    print("  O(modes^2) in mode-to-mode work, so the MEASURED ratios SHOULD EXCEED these.")
    print("  A reader who mistakes the model for an expectation will misread that as a defect.")


def print_extras(legs):
    """SeisSol-only cross-checks that cost nothing and audit the headline."""
    rows = []
    for p in ORDERS:
        leg = legs["seissol_p%d" % p]
        if leg.log is None:
            continue
        log = leg.log
        hw, nz = log.get("hw_flop"), log.get("nz_flop")
        rows.append((leg.label,
                     log.get("load_imbalance_pct"),
                     hw, nz,
                     (100.0 * nz / hw) if (hw and nz) else None,
                     (hw / leg.m["t_loop"] / leg.m["cores_active"] / 1e9)
                     if (hw and leg.m.get("t_loop") and leg.m.get("cores_active"))
                     else None,
                     log.get("per_element")))
    if not rows:
        return
    print("")
    print(rule("="))
    print("6. SeisSol-only cross-checks  (free; they audit the headline, not replace it)")
    print(rule("="))
    print("  %-22s %10s %14s %14s %9s %14s"
          % ("leg", "imbal %", "HW-FLOP", "NZ-FLOP", "NZ/HW %", "GFLOP/s/core"))
    for lbl, imb, hw, nz, pct, gfs, _pe in rows:
        print("  %-22s %10s %14s %14s %9s %14s"
              % (lbl, fmt(imb, "%.2f"), fmt(hw, "%.4e"), fmt(nz, "%.4e"),
                 fmt(pct, "%.1f"), fmt(gfs, "%.2f")))
    print("")
    print("  GFLOP/s/core recomputed against 'Simulation time (compute)' -- SeisSol's own")
    print("  printPerformanceSummary divides by compute+blocking-IO (Simulator.cpp:148,163).")
    print("  Low NZ/HW at o2 = vector-lane padding of the 4-mode kernels; it is the reason")
    print("  the p1 cross-code ratio is the least reliable of the three.")
    print("")
    print("  LoopStatistics '( per element )' regression slopes (seconds per element-update")
    print("  per RANK-call, LoopStatistics.cpp:226-236).  These are a DIFFERENT denominator")
    print("  (per rank-team of OMP threads, kernel-only) -- use them only as an internal")
    print("  consistency check on SeisSol's own T_loop, never against the MFEM number:")
    for lbl, _i, _h, _n, _p, _g, pe in rows:
        if pe:
            print("    %-22s %s" % (lbl, ", ".join("%s=%.3e" % (k, v)
                                                   for k, v in sorted(pe.items()))))
        else:
            print("    %-22s (no regression lines parsed)" % lbl)


def print_uncomputed(legs):
    print("")
    print(rule("="))
    print("7. WHAT THIS RUN COULD **NOT** COMPUTE   (printed, never silently omitted)")
    print(rule("="))
    n_required = 0
    any_note = False
    for code in ("mfem", "seissol"):
        for p in ORDERS:
            leg = legs["%s_p%d" % (code, p)]
            if not leg.missing:
                continue
            any_note = True
            print("")
            print("  %s" % leg.label)
            for q, why, required in leg.missing:
                if required:
                    n_required += 1
                print("    - [%s] %-28s : %s"
                      % ("BLOCKING " if required else "disclosed", q, why))
    if not any_note:
        print("  (nothing -- every quantity required by R6 was obtained for all six legs)")
    print("")
    print("  [BLOCKING ] = an R6 quantity that is genuinely absent; the affected cell is")
    print("                blank above and that leg is NOT part of the baseline.")
    print("  [disclosed] = a quantity the code structurally cannot emit; it was substituted")
    print("                from the pre-registered contract AND gated, so the leg stands.")
    return n_required


def print_caveat():
    print("")
    print(rule("="))
    print("CAVEAT (one line, per the recon -- reproduce it verbatim wherever this table is quoted):")
    print("  This is a COST baseline at NOMINALLY-matched order on ONE mesh, ONE node count,")
    print("  ONE problem, GTS only, one run per leg -- it measures NO accuracy, so it cannot")
    print("  claim 'time to solution at matched accuracy', any Pareto/'p3 is worth it' result,")
    print("  any scaling behaviour, and nothing about LTS or any other mesh.")
    print(rule("="))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv):
    ap = argparse.ArgumentParser(
        description="TPV104 200 m GTS order baseline post-processor "
                    "(MFEM p1/p2/p3 vs SeisSol o2/o3/o4).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="directory holding the six legs (scanned recursively)")
    ap.add_argument("--tfinal", type=float, default=DEFAULT_TFINAL,
                    help="simulated window in seconds (default %(default)s)")
    ap.add_argument("--elements", type=int, default=PINNED_ELEMENTS,
                    help="mesh element count (default %(default)s)")
    ap.add_argument("--dt-ref-p3", type=float, default=DT_REF_P3,
                    help="pinned p3/o4 dt in seconds (default %(default)s)")
    ap.add_argument("--mfem-region", default=DEFAULT_REGION,
                    help="Caliper region for the MFEM time loop (default %(default)s)")
    ap.add_argument("--mfem-cores", type=int, default=None,
                    help="override MFEM cores per leg (default: the 'ranks:' banner)")
    ap.add_argument("--seissol-cores", type=int, default=None,
                    help="override SeisSol ALLOCATED cores per leg")
    ap.add_argument("--seissol-cpus-per-rank", type=int,
                    default=DEFAULT_SEISSOL_CPUS_PER_RANK,
                    help="cores allocated per SeisSol rank (default %(default)s)")
    ap.add_argument("--dt-tol", type=float, default=1e-5,
                    help="relative dt agreement tolerance (default %(default)s)")
    ap.add_argument("--json", default=None, help="also write results as JSON")
    ap.add_argument("--strict", action="store_true",
                    help="also return rc 3 when only WARN gates fired "
                         "(FAIL gates already return 3)")
    args = ap.parse_args(argv[1:])

    _set_refs(args.dt_ref_p3, args.mfem_region)

    if not os.path.isdir(args.run_dir):
        sys.stderr.write("ERROR: not a directory: %s\n" % args.run_dir)
        return 2
    if args.tfinal <= 0:
        sys.stderr.write("ERROR: --tfinal must be > 0\n")
        return 2

    legs, stats = build_legs(args.run_dir, args.mfem_region)

    print(rule("="))
    print("TPV104 200 m -- GTS ORDER SPEED BASELINE  (MFEM p1/p2/p3 vs SeisSol o2/o3/o4)")
    print(rule("="))
    print("  Fairness contract enforced here: R1 same mesh/EndTime/CFL/GTS | R2 accuracy-")
    print("  matched p_k<->o_(k+1) | R3 per-order dt, both s/sim-s AND per-element-update")
    print("  reported | R4 IO off on both sides | R5 each code in ITS validated rank shape")
    print("  | R6 every leg reports wall, steps, dt, elements, cores.")

    print_contract(args)
    missing_legs = print_inventory(args.run_dir, legs, stats)

    for leg in legs.values():
        derive(leg, args)

    ss_elems = None
    for p in ORDERS:
        e = legs["seissol_p%d" % p].log
        if e and e.get("elements"):
            ss_elems = e["elements"]
            break
    for leg in legs.values():
        check_gates(leg, args, ss_elems)

    # cross-code dt equality (R1/R2): the single cheapest fairness check there is
    print("")
    print(rule("="))
    print("2b. CROSS-CODE dt / STEP EQUALITY  (the check that validates the whole comparison)")
    print(rule("="))
    dt_fail = 0
    for p in ORDERS:
        a = legs["mfem_p%d" % p].m
        b = legs["seissol_p%d" % p].m
        da, db = a.get("dt"), b.get("dt")
        if da is None or db is None:
            print("  p%d <-> o%d : dt UNKNOWN on %s side -- equality UNVERIFIED"
                  % (p, p + 1, "the MFEM" if da is None else "the SeisSol"))
            dt_fail += 1
            continue
        rel = abs(da - db) / max(da, db)
        ok = rel <= args.dt_tol
        if not ok:
            dt_fail += 1
        print("  p%d <-> o%d : MFEM %.6e s   SeisSol %.6e s   rel dev %.2e   [%s]"
              % (p, p + 1, da, db, rel, "PASS" if ok else "FAIL"))
        sa, sb = a.get("steps"), b.get("steps")
        if sa and sb:
            print("              steps: MFEM %d, SeisSol %d %s"
                  % (sa, sb, "" if sa == sb else "  <-- MISMATCH, the s/sim-s "
                                                 "ratio is no longer a pure per-step ratio"))

    n_fail = print_gates(legs)
    print_baseline(legs, args)
    print_cross_code(legs, args)
    print_within_code(legs, args)
    print_extras(legs)
    n_blocking = print_uncomputed(legs)
    print_caveat()

    print("")
    print("SUMMARY: %d/6 legs found, %d blocking gap(s), %d gate FAIL(s), "
          "%d cross-code dt check FAIL(s)."
          % (sum(1 for l in legs.values() if l.log is not None),
             n_blocking, n_fail, dt_fail))

    if args.json:
        payload = {
            "run_dir": os.path.abspath(args.run_dir),
            "tfinal": args.tfinal,
            "elements": args.elements,
            "dt_ref_p3": args.dt_ref_p3,
            "contract": {"p%d" % p: {"dt": dt_pinned(p),
                                     "steps": steps_pinned(p, args.tfinal),
                                     "modes": MODES[p],
                                     "ader_order": p + 1,
                                     "seissol_order": p + 1}
                         for p in ORDERS},
            "legs": {},
            "gate_failures": n_fail,
            "dt_check_failures": dt_fail,
        }
        for key, leg in sorted(legs.items()):
            payload["legs"][key] = {
                "label": leg.label,
                "log": leg.log["path"] if leg.log else None,
                "caliper": leg.cali["path"] if leg.cali else None,
                "deck": leg.deck["path"] if leg.deck else None,
                "metrics": {k: v for k, v in leg.m.items()
                            if isinstance(v, (int, float, str)) or v is None},
                "gates": [{"name": n, "status": s, "detail": d}
                          for n, s, d in leg.gates],
                "missing": [{"quantity": q, "why": w, "required": r}
                            for q, w, r in leg.missing],
            }
        try:
            with open(args.json, "w") as fh:
                json.dump(payload, fh, indent=2, sort_keys=True)
            print("JSON written to %s" % args.json)
        except (OSError, IOError) as ex:
            sys.stderr.write("WARNING: could not write JSON (%s)\n" % ex)

    if missing_legs or n_blocking:
        return 1
    if n_fail or dt_fail:
        return 3
    return 3 if (args.strict and any(s == "WARN" for l in legs.values()
                                     for _, s, _ in l.gates)) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
