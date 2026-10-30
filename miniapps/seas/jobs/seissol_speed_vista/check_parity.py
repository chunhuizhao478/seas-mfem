#!/usr/bin/env python3
"""Physics parity check: MFEM Nitsche deck  vs  SeisSol deck (this folder).

Usage:
  python3 check_parity.py --problem tpv205|tpv104|tpv102 --repo <repo root>

MFEM side: the toml the MFEM speed runs use (mfem_decks/<prob>_nitsche_sst.toml,
checked equal to its source deck except output.paraview_enabled), read with
tomllib, plus the driver's code-level rules, re-implemented here from the
C++ with file:line comments (paths relative to miniapps/seas/):
  prestress     drivers/nitsche_dyn_driver.cpp:1567-1645 (global tensor,
                [[stress_region]] closed boxes, last match wins)
  pin frame     nitsche_dyn_driver.cpp:1647-1654 (pn = ref_normal,
                ps = up x pn / |.|, pd = ps x pn)
  nucleation    nitsche_dyn_driver.cpp:1656-1684 (dsigma = -dtau F G
                (s n + n s), r measured with (pd, ps)); F = CompactBellFactor
                dynamic/spatial_nucleation.cpp:258-265; G = SmoothStep :20-29
  RS fields     spatial/code/spatial_friction.cpp:2146-2205 (seed defaults,
                boxcar_taper param = outer + (inner - outer) B),
                BoxcarTaperFactor :111-124, SCECBoxcar
                spatial/code/spatial_friction.hpp:336-342
  LSW fields    spatial/code/spatial_friction.cpp:1999-2010 (defaults)
  psi_0         friction/dieterich_ruina.hpp:507-529 (InitialStatePsi, eta 0),
                seeded at nitsche_dyn_driver.cpp:1812-1835
  station 9th   nitsche_dyn_driver.cpp:2005-2017
SeisSol side: parameters.par (Fortran namelist, keys lower-cased as fty does),
<prob>_fault.yaml / _material.yaml (PyYAML; easi tags kept as data; LuaMap
bodies translated to Python by the small translator below and evaluated),
and SeisSol's own formulas re-implemented with file:line (commit a8e7f04f):
  psi_0 FL=103  DynamicRupture/Initializer/RateAndStateInitializer.cpp:154-180
  theta_0 FL=3  RateAndStateInitializer.cpp:104-132 (psi = f0 + b ln(V0 theta/L))
  logsinh       DynamicRupture/FrictionLaws/RateAndStateCommon.h:135-139
  smoothStep    Numerical/GaussianNucleationFunction.h:22-39
  strike/dip    DynamicRupture/Misc.cpp:19-34; plus side Geometry/MeshReader.cpp:139-163
  Ts0/Td0/Pn0   DynamicRupture/Output/ReceiverBasedOutput.cpp:308-325
Fields are compared at ~200 fault points, relative tolerance 1e-9.
Exit status: 0 = every compared item matches, 1 = a mismatch, 2 = input error.
Items that are not compared (no SeisSol counterpart) are printed as OPEN and
do not change the exit status; NOTES.md lists them.
"""
import argparse, math, os, re, sys
import tomllib

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML missing: python3 -m pip install pyyaml")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REL = 1e-9

PROBLEMS = {
    "tpv205": dict(deck="miniapps/seas/jobs/seissol_speed_vista/mfem_decks/tpv205_nitsche_sst.toml",
                   source_deck="miniapps/seas/tpv205/configs/tpv205_nitsche.toml",
                   stations="miniapps/seas/tpv205/stations/tpv205_fault_stations.csv",
                   mesh_geo="miniapps/seas/tpv205/mesh/tpv2053d_200m.geo",
                   meshfile="tpv2053d_200m.puml.h5",
                   fault_box=(-15000.0, 15000.0, -15000.0, 0.0)),
    "tpv104": dict(deck="miniapps/seas/jobs/seissol_speed_vista/mfem_decks/tpv104_nitsche_sst.toml",
                   source_deck="miniapps/seas/tpv104/configs/tpv104_nitsche.toml",
                   stations="miniapps/seas/tpv104/stations/tpv104_fault_stations.csv",
                   mesh_geo="miniapps/seas/tpv104/mesh/tpv104_200m.geo",
                   meshfile="tpv104_200m.puml.h5",
                   fault_box=(-18000.0, 18000.0, -18000.0, 0.0)),
    "tpv102": dict(deck="miniapps/seas/jobs/seissol_speed_vista/mfem_decks/tpv102_nitsche_sst.toml",
                   source_deck="miniapps/seas/tpv102/configs/tpv102_nitsche_pc_p2_kv0p1.toml",
                   stations="miniapps/seas/tpv102/stations/tpv102_fault_stations.csv",
                   mesh_geo="miniapps/seas/tpv104/mesh/tpv104_200m.geo",
                   meshfile="tpv104_200m.puml.h5",
                   fault_box=(-18000.0, 18000.0, -18000.0, 0.0)),
}

# MFEM off-fault (free-surface) stations: the Nitsche driver writes none; the
# native setups list these (x, fault-normal y, z = 0), in this order.
SURFACE_STATIONS = {
    # dynamic/tpv205_setup.hpp:195-210: intentionally empty
    "tpv205": ([], "dynamic/tpv205_setup.hpp:207-210"),
    # dynamic/tpv104_setup.hpp:209-219
    "tpv104": ([(0.0, 9e3, 0.0), (0.0, -9e3, 0.0), (12e3, 6e3, 0.0), (12e3, -6e3, 0.0),
                (-12e3, 6e3, 0.0), (-12e3, -6e3, 0.0)], "dynamic/tpv104_setup.hpp:209-219"),
    # dynamic/tpv102_setup.hpp:178-189
    "tpv102": ([(0.0, 9e3, 0.0), (0.0, -9e3, 0.0), (12e3, 6e3, 0.0), (12e3, -6e3, 0.0),
                (-12e3, 6e3, 0.0), (-12e3, -6e3, 0.0)], "dynamic/tpv102_setup.hpp:178-189"),
}

# SeisSol default face map, src/Initializer/FaceMap.cpp:98-117
def seissol_face_type(tag):
    if tag in (0, 6):
        return "regular"
    if tag == 1:
        return "free"
    if tag == 2:
        return "free_gravity"
    if tag == 3 or tag >= 65:
        return "fault"
    if tag == 4:
        return "dirichlet"
    if tag == 5:
        return "absorbing"
    if tag == 7:
        return "analytical"
    return "unmapped"


# --------------------------------------------------------------------- report
class Report:
    def __init__(self):
        self.rows, self.bad, self.open = [], 0, []

    def add(self, name, mfem, mref, ss, sref, ok):
        self.rows.append((name, mfem, mref, ss, sref, "OK" if ok else "MISMATCH"))
        if not ok:
            self.bad += 1

    def scalar(self, name, mfem, mref, ss, sref):
        self.add(name, fmt(mfem), mref, fmt(ss), sref, close(mfem, ss))

    def info(self, name, mfem, mref, ss, sref):
        self.open.append((name, mfem, mref, ss, sref))

    def show(self):
        """Side by side: quantity | MFEM value | SeisSol value | status, then
        the two sources on an indented second line."""
        w0, w1 = 36, 30
        print(f"{'quantity':<{w0}} | {'MFEM':<{w1}} | {'SeisSol':<{w1}} | status")
        print("-" * (w0 + 2 * w1 + 16))
        for name, mv, mref, sv, sref, st in self.rows:
            print(f"{name:<{w0}} | {str(mv):<{w1}} | {str(sv):<{w1}} | {st}")
            print(f"{'':<{w0}}   src MFEM: {mref}  |  SeisSol: {sref}")
        if self.open:
            print("\nOPEN (no SeisSol counterpart in this deck; not compared, see NOTES.md):")
            for r in self.open:
                print(f"  {r[0]}: MFEM {r[1]} ({r[2]}) | SeisSol {r[3]} ({r[4]})")


def fmt(v):
    if isinstance(v, float):
        return f"{v:.10g}"
    return str(v)


def close(a, b, rel=REL):
    if isinstance(a, (int, float, np.floating)) and isinstance(b, (int, float, np.floating)):
        a, b = float(a), float(b)
        if a == b:
            return True
        if not (math.isfinite(a) and math.isfinite(b)):
            return False
        return abs(a - b) <= rel * max(abs(a), abs(b))
    return a == b


# ------------------------------------------------------------------ toml side
def toml_line(path, section, key):
    """1-based line of `key =` inside [section] (or [[section]] occurrence k)."""
    sec, occ = (section, 1) if not isinstance(section, tuple) else section
    cur, seen = None, 0
    for i, line in enumerate(open(path), 1):
        s = line.split("#")[0].strip()
        m = re.match(r"^\[\[?([^\]]+)\]\]?$", s)
        if m:
            cur = m.group(1).strip()
            if cur == sec:
                seen += 1
            continue
        if cur == sec and seen == occ and re.match(rf"^{re.escape(key)}\s*=", s):
            return i
    return None


def seconds(v):
    if isinstance(v, (int, float)):
        return float(v)
    m = re.fullmatch(r"\s*([-+0-9.eE]+)\s*s\s*", str(v))
    if not m:
        raise ValueError(f"cannot read time {v!r}")
    return float(m.group(1))


class MfemDeck:
    def __init__(self, path):
        self.path = path
        with open(path, "rb") as f:
            self.t = tomllib.load(f)

    def get(self, section, key):
        node = self.t
        for p in section.split("."):
            node = node[p]
        return node[key]

    def ref(self, section, key, occ=1):
        ln = toml_line(self.path, (section, occ), key)
        return f"{os.path.basename(self.path)}:{ln}"


def mfem_frame(deck):
    """nitsche_dyn_driver.cpp:1647-1654."""
    pn = np.array(deck.get("fault_geometry", "ref_normal"), float)
    pu = np.array(deck.get("fault_geometry", "up"), float)
    pn = pn / np.linalg.norm(pn)
    ps = np.cross(pu, pn)
    ps = ps / np.linalg.norm(ps)
    pd = np.cross(ps, pn)
    return pn, ps, pd


def mfem_prestress(deck, x):
    """Global tensor at x (3x3, tension positive); driver :1567-1645."""
    st = deck.t["stress"]
    assert st["kind"] == "constant_tensor"

    def tens(d):
        xx, yy, zz = d["sigma_xx_pa"], d["sigma_yy_pa"], d["sigma_zz_pa"]
        xy, yz, xz = d["sigma_xy_pa"], d["sigma_yz_pa"], d["sigma_xz_pa"]
        return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], float)

    s = tens(st)
    for r in deck.t.get("stress_region", []):
        lo = [r.get(f"{c}_min_m", -math.inf) for c in "xyz"]
        hi = [r.get(f"{c}_max_m", math.inf) for c in "xyz"]
        if all(lo[i] <= x[i] <= hi[i] for i in range(3)):
            s = tens(r)
    return s


def compact_bell(r, R):
    """dynamic/spatial_nucleation.cpp:258-265."""
    if r >= R:
        return 0.0
    r2, R2 = r * r, R * R
    return math.exp(r2 / (r2 - R2))


def mfem_smoothstep(t, t0):
    """dynamic/spatial_nucleation.cpp:20-29."""
    if t <= 0.0:
        return 0.0
    if t < t0:
        tau = t - t0
        return math.exp(tau * tau / (t * (t - 2.0 * t0)))
    return 1.0


def mfem_nucleation(deck, x):
    """dsigma (3x3) per unit G, driver :1656-1684; None without [nucleation]."""
    if "nucleation" not in deck.t:
        return None
    nk = deck.t["nucleation"]["kind"]
    assert nk == "gradual_overstress_compact_circular", nk
    sp = deck.t["nucleation"][nk]
    pn, ps, pd = mfem_frame(deck)
    c = np.array([sp["center_x_m"], sp["center_y_m"], sp["center_z_m"]], float)
    e = np.asarray(x, float) - c
    ed, es = float(e @ pd), float(e @ ps)
    amp = -sp["delta_tau_pa"] * compact_bell(math.sqrt(ed * ed + es * es), sp["radius_m"])
    return amp * (np.outer(ps, pn) + np.outer(pn, ps))


def scec_boxcar(offset, half, trans):
    """spatial/code/spatial_friction.hpp:336-342."""
    ax = abs(offset)
    if ax <= half:
        return 1.0
    if trans <= 0.0 or ax >= half + trans:
        return 0.0
    return 0.5 * (1.0 + math.tanh(trans / (ax - half - trans) + trans / (ax - half)))


def boxcar_factor(rule, x):
    """spatial_friction.cpp:111-124 (non-finite half = untapered axis)."""
    f = 1.0
    for i, c in enumerate("xyz"):
        half = rule.get(f"boxcar_half_{c}_m", math.inf)
        if not math.isfinite(half):
            continue
        f *= scec_boxcar(x[i] - rule.get(f"boxcar_center_{c}_m", 0.0), half,
                         rule.get(f"boxcar_trans_{c}_m", 0.0))
    return f


def mfem_rs_fields(deck, x):
    """a, V_w, Dc, V_init per point; spatial_friction.cpp:2146-2205."""
    rs = deck.t["friction"]["rate_state"]
    a, V_w = rs["a_default"], rs.get("V_w_default", math.nan)
    for r in rs.get("spatial", []):
        if r["kind"] == "boxcar_taper":
            f = boxcar_factor(r, x)
            if "a_inner" in r:
                a = r["a_outer"] + (r["a_inner"] - r["a_outer"]) * f
            if "V_w_inner" in r:
                V_w = r["V_w_outer"] + (r["V_w_inner"] - r["V_w_outer"]) * f
        else:
            raise NotImplementedError(f"rule kind {r['kind']} not reimplemented")
    return dict(a=a, V_w=V_w, Dc=rs["Dc_default"], V_init=rs["V_init_default"])


def mfem_psi0(deck, x):
    """InitialStatePsi (dieterich_ruina.hpp:507-529, eta 0) at the t = 0 traction."""
    pn, ps, pd = mfem_frame(deck)
    S = mfem_prestress(deck, x)            # G(0) = 0: no nucleation at t = 0
    T = S @ pn
    sig = -float(T @ pn)
    X = math.hypot(float(T @ pd), float(T @ ps))
    rs = deck.t["friction"]["rate_state"]
    f = mfem_rs_fields(deck, x)
    fa = (X / sig) / f["a"]
    sinh_val = 0.5 * math.exp(fa) if fa > 700.0 else math.sinh(fa)
    return f["a"] * math.log((2.0 * rs["V_0_default"] / f["V_init"]) * sinh_val)


# --------------------------------------------------------------- SeisSol side
def read_par(path):
    """Fortran namelist -> {section: {key: (value string, line)}} (keys lower-case)."""
    out, sec = {}, None
    for i, raw in enumerate(open(path), 1):
        line = raw.split("!")[0].strip()
        if not line:
            continue
        if line.startswith("&"):
            sec = line[1:].strip().lower()
            out.setdefault(sec, {})
            continue
        if line == "/":
            sec = None
            continue
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$", line)
        if not m or sec is None:
            raise ValueError(f"{path}:{i}: cannot parse {raw!r}")
        v = m.group(2).strip().strip("'").strip('"')
        out[sec][m.group(1).lower()] = (v, i)
    return out


def par_num(v):
    return float(v.replace("d", "e").replace("D", "e"))


class _Loader(yaml.SafeLoader):
    pass


def _mapping(loader, node, deep=True):
    out = {}
    for kn, vn in node.value:
        k = loader.construct_object(kn, deep=True)
        if isinstance(k, list):
            k = tuple(k)
        out[k] = loader.construct_object(vn, deep=True)
    return out


def _tagged(loader, suffix, node):
    if isinstance(node, yaml.MappingNode):
        body = _mapping(loader, node)
    elif isinstance(node, yaml.SequenceNode):
        body = loader.construct_sequence(node, deep=True)
    else:
        body = loader.construct_scalar(node)
    return {"__tag__": "!" + suffix, "body": body}


_Loader.add_multi_constructor("!", _tagged)
_Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def load_easi(path):
    with open(path) as f:
        return yaml.load(f, Loader=_Loader)


class LuaMath:
    abs = staticmethod(abs)
    max = staticmethod(max)
    min = staticmethod(min)
    exp = staticmethod(lambda v: np.exp(v))
    sqrt = staticmethod(lambda v: np.sqrt(v))
    tanh = staticmethod(lambda v: np.tanh(v))
    log = staticmethod(lambda v: np.log(v))
    pi = math.pi


def lua_to_python(src):
    """Translate the Lua subset the decks use (function f (x), assignments,
    if/elseif/else/end, return { k = e, ... }) to Python source.  Anything
    else raises, so a change in the yaml cannot be evaluated silently wrong."""
    def expr(e):
        e = e.replace("~=", "!=").replace("^", "**")
        e = re.sub(r"\bmath\.", "_math.", e)
        if re.search(r"\.\.|#|\bnil\b|\blocal\b|:", e):
            raise ValueError(f"unsupported Lua expression {e!r}")
        return e

    out, ind, ret = [], 0, None
    for raw in src.splitlines():
        s = raw.split("--")[0].strip()
        if not s:
            continue
        if ret is not None:
            if s.startswith("}"):
                out.append("    " * ind + "return {" + ", ".join(ret) + "}")
                ret = None
                continue
            m = re.match(r"^([A-Za-z_]\w*)\s*=\s*(.+?),?$", s)
            if not m:
                raise ValueError(f"bad return field {s!r}")
            ret.append(f"{m.group(1)!r}: {expr(m.group(2))}")
            continue
        if m := re.match(r"^function\s+f\s*\(\s*x\s*\)$", s):
            out.append("def f(x):"); ind = 1; continue
        if m := re.match(r"^if\s+(.+)\s+then$", s):
            out.append("    " * ind + f"if {expr(m.group(1))}:"); ind += 1; continue
        if m := re.match(r"^elseif\s+(.+)\s+then$", s):
            out.append("    " * (ind - 1) + f"elif {expr(m.group(1))}:"); continue
        if s == "else":
            out.append("    " * (ind - 1) + "else:"); continue
        if s == "end":
            ind -= 1
            if ind < 0:
                raise ValueError("unbalanced end")
            continue
        if s == "return {":
            ret = []; continue
        if m := re.match(r"^([A-Za-z_]\w*)\s*=\s*(.+)$", s):
            out.append("    " * ind + f"{m.group(1)} = {expr(m.group(2))}"); continue
        raise ValueError(f"unsupported Lua line {s!r}")
    if ind != 0 or ret is not None:
        raise ValueError("unbalanced Lua body")
    return "\n".join(out)


class EasiModel:
    """Evaluate a top-level !Switch of !ConstantMap / !LuaMap at a point."""
    def __init__(self, path):
        self.path = path
        doc = load_easi(path)
        self.branches = []
        if isinstance(doc, dict) and doc.get("__tag__") == "!Switch":
            items = doc["body"].items()
        elif isinstance(doc, dict) and doc.get("__tag__") == "!ConstantMap":
            items = [(tuple(doc["body"]["map"].keys()), doc)]
        else:
            raise ValueError(f"{path}: top level must be !Switch or !ConstantMap")
        for keys, node in items:
            keys = keys if isinstance(keys, tuple) else (keys,)
            tag = node["__tag__"]
            if tag == "!ConstantMap":
                mp = {k: float(v) for k, v in node["body"]["map"].items()}
                self.branches.append((keys, "const", mp))
            elif tag == "!LuaMap":
                py = lua_to_python(node["body"]["function"])
                ns = {"_math": LuaMath}
                exec(py, ns)
                self.branches.append((keys, "lua", (ns["f"], tuple(node["body"]["returns"]))))
            else:
                raise ValueError(f"{path}: easi tag {tag} not handled")
        self.provides = {k for keys, _, _ in self.branches for k in keys}

    def __call__(self, name, x):
        for keys, kind, data in self.branches:
            if name in keys:
                if kind == "const":
                    return data[name]
                f, rets = data
                with np.errstate(all="ignore"):
                    r = f({"x": np.float64(x[0]), "y": np.float64(x[1]), "z": np.float64(x[2])})
                return float(r[name])
        raise KeyError(f"{self.path}: no branch provides {name}")


def seissol_smoothstep(t, t0):
    """Numerical/GaussianNucleationFunction.h:22-31."""
    if t <= 0:
        return 0.0
    if t < t0:
        tau = t - t0
        return math.exp(tau * tau / (t * (t - 2.0 * t0)))
    return 1.0


def logsinh(x, c):
    """RateAndStateCommon.h:135-139."""
    sign = 1.0 if c >= 0 else -1.0
    ac = abs(c)
    return ac + math.log(x / 2 * -sign * math.expm1(-2 * ac))


def seissol_frame(par):
    """Plus side and fault normal of the y = 0 fault (MeshReader.cpp:139-163),
    strike/dip (Misc.cpp:19-34).  Returns n, strike, dip."""
    dr = par["dynamicrupture"]
    ref = np.array([par_num(dr.get(k, ("0", 0))[0]) for k in ("xref", "yref", "zref")])
    method = int(dr.get("refpointmethod", ("0", 0))[0])
    # outward normal of the y<0 cell is +y; of the y>0 cell is -y
    if method == 0:   # plus cell = cell on the same side as the reference point
        n = np.array([0.0, 1.0, 0.0]) if ref[1] < 0 else np.array([0.0, -1.0, 0.0])
    else:             # plus cell = cell whose outward normal has dot(ref, n) > 0
        n = np.array([0.0, 1.0, 0.0]) if ref[1] > 0 else np.array([0.0, -1.0, 0.0])
    inv = 1.0 / math.sqrt(n[0] ** 2 + n[1] ** 2)
    s = np.array([n[1] * inv, -n[0] * inv, 0.0])
    d = np.array([-s[1] * n[2], s[0] * n[2], s[1] * n[0] - s[0] * n[1]])
    d = d / np.linalg.norm(d)
    return n, s, d


# ------------------------------------------------------------------- points
def sample_points(box, stations, extra):
    x0, x1, z0, z1 = box
    rng = np.random.default_rng(20261030)
    pts = [(float(x), 0.0, float(z)) for x, z in zip(rng.uniform(x0, x1, 200), rng.uniform(z0, z1, 200))]
    for x in np.linspace(x0, x1, 7):
        for z in np.linspace(z0, z1, 7):
            pts.append((float(x), 0.0, float(z)))
    pts += [tuple(map(float, s)) for s in stations]
    pts += [p for p in extra if x0 <= p[0] <= x1 and z0 <= p[2] <= z1]
    return pts


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", required=True, choices=sorted(PROBLEMS))
    ap.add_argument("--repo", required=True)
    a = ap.parse_args()
    P = PROBLEMS[a.problem]
    repo = os.path.abspath(a.repo)
    deck = MfemDeck(os.path.join(repo, P["deck"]))
    ddir = os.path.join(HERE, "decks", a.problem)
    par_path = os.path.join(ddir, "parameters.par")
    par = read_par(par_path)
    pr = lambda sec, key: f"parameters.par:{par[sec][key][1]}"
    fault = EasiModel(os.path.join(ddir, f"{a.problem}_fault.yaml"))
    mat = EasiModel(os.path.join(ddir, f"{a.problem}_material.yaml"))
    fy = f"{a.problem}_fault.yaml"
    R = Report()
    print(f"problem {a.problem}\n  MFEM deck   {P['deck']}\n  SeisSol dir {os.path.relpath(ddir, repo)}\n")

    # ---- the speed-test MFEM deck (what run_pair_vista.sh runs) is its source
    # deck with only the output switch paraview_enabled changed
    with open(os.path.join(repo, P["source_deck"]), "rb") as f:
        src_t = tomllib.load(f)
    strip = lambda t: {k: ({kk: vv for kk, vv in v.items() if kk != "paraview_enabled"}
                           if isinstance(v, dict) else v) for k, v in t.items()}
    R.add("MFEM speed deck = source deck (but paraview_enabled)", "same parsed toml",
          f"{P['deck']} vs {P['source_deck']}", "-", "-", strip(src_t) == strip(deck.t))

    # ---- law / FL
    law = deck.get("meta", "law")
    evo = deck.t.get("friction", {}).get("rate_state", {}).get("state_evolution")
    want_fl = {"slip_weakening": 16}.get(law) if law == "slip_weakening" else \
        {"aging_law": 3, "slip_law_strong_rate_weakening": 103}.get(evo)
    fl = int(par["dynamicrupture"]["fl"][0])
    R.add("friction law", f"{law}/{evo}" if evo else law,
          deck.ref("meta", "law") + (", " + deck.ref("friction.rate_state", "state_evolution") if evo else ""),
          f"FL={fl}", pr("dynamicrupture", "fl") + " (Misc.h:33-47)", fl == want_fl)

    # ---- material
    mc = deck.t["material_constant_fallback"]
    for k in ("rho", "mu", "lambda"):
        R.scalar(f"material {k}", float(mc[k]), deck.ref("material_constant_fallback", k),
                 mat(k, (0, 0, 0)), f"{a.problem}_material.yaml")

    # ---- time
    R.scalar("end time [s]", seconds(deck.get("time", "tfinal")), deck.ref("time", "tfinal"),
             par_num(par["abortcriteria"]["endtime"][0]), pr("abortcriteria", "endtime"))

    # ---- boundary tags
    b = deck.t["boundary"]
    remap = {101: 1, 103: 3, 105: 5} if a.problem == "tpv205" else {}
    for role, tags in (("fault", [b["fault_attr"]]), ("free", b["natural_attrs"]),
                       ("absorbing", b["absorbing_attrs"])):
        for t in tags:
            tt = remap.get(t, t)
            ok = seissol_face_type(tt) == role
            src = "FaceMap.cpp:98-117" + (f" after retag_msh2.py {t}->{tt}" if remap else "")
            R.add(f"tag {t} ({role})", role, deck.ref("boundary", {"fault": "fault_attr", "free":
                  "natural_attrs", "absorbing": "absorbing_attrs"}[role]),
                  f"{tt} -> {seissol_face_type(tt)}", src, ok)
    R.add("mesh file", os.path.basename(P["mesh_geo"]), P["mesh_geo"],
          par["meshnml"]["meshfile"][0], pr("meshnml", "meshfile"),
          par["meshnml"]["meshfile"][0] == P["meshfile"])

    # ---- stations
    st_csv = os.path.join(repo, P["stations"])
    stations = []
    for line in open(st_csv):
        if line.startswith("#") or line.startswith("name"):
            continue
        n_, x_, y_, z_ = line.strip().split(",")
        stations.append((float(x_), float(y_), float(z_)))
    fr = [tuple(float(v) for v in l.split()) for l in open(os.path.join(ddir, f"{a.problem}_faultreceivers.dat")) if l.strip()]
    R.add("on-fault stations (order, xyz)", f"{len(stations)} pts", os.path.basename(st_csv),
          f"{len(fr)} pts", f"{a.problem}_faultreceivers.dat", fr == stations)
    surf, sref = SURFACE_STATIONS[a.problem]
    rc = [tuple(float(v) for v in l.split()) for l in open(os.path.join(ddir, f"{a.problem}_receivers.dat")) if l.strip()]
    R.add("off-fault receivers (order, xyz)", f"{len(surf)} pts", sref, f"{len(rc)} pts",
          f"{a.problem}_receivers.dat", rc == surf)
    out = par["output"]
    want_rc = "1" if surf else "0"     # TPV205: no off-fault stations in MFEM, so the receivers are off
    R.add("off-fault receiver output", f"{'on' if surf else 'off'} (pickdt 0.005, task)", "task spec",
          f"ReceiverOutput={out['receiveroutput'][0]} pickdt={out['pickdt'][0]}",
          pr("output", "receiveroutput"), out["receiveroutput"][0] == want_rc and par_num(out["pickdt"][0]) == 0.005)
    R.add("on-fault output", "pickpoints only", "task spec",
          f"OutputPointType={par['dynamicrupture']['outputpointtype'][0]}",
          pr("dynamicrupture", "outputpointtype") + " (DRParameters.h:22-27)",
          par["dynamicrupture"]["outputpointtype"][0] == "3")

    # ---- points: random + grid (sample_points), dense lines through the
    # hypocentre and the boxcar transition, and points 1 m / 50 m on either
    # side of every edge the fields have (patch boxes, boxcar plateau and
    # transition ends, bell radius)
    extra = []
    edges_x = (0.0, 1500.0, 6000.0, 9000.0, 15000.0, 18000.0)
    edges_z = (0.0, -6000.0, -9000.0, -15000.0, -18000.0, -7500.0)
    for ex in edges_x:
        for sx in (1.0, -1.0):
            for dx in (0.0, 1.0, -1.0, 50.0, -50.0):
                for ez in edges_z:
                    for dz in (0.0, 1.0, -1.0, 50.0, -50.0):
                        extra.append((sx * ex + dx, 0.0, ez + dz))
    x0, x1, z0, z1 = P["fault_box"]
    for xv in np.arange(x0, x1 + 1.0, 250.0):
        extra += [(float(xv), 0.0, -7500.0), (float(xv), 0.0, -16500.0), (float(xv), 0.0, -100.0)]
    for zv in np.arange(z0, z1 + 1.0, 250.0):
        extra += [(0.0, 0.0, float(zv)), (16500.0, 0.0, float(zv)), (-7500.0, 0.0, float(zv))]
    for r in (0.0, 1000.0, 2999.0, 2999.999, 3000.0, 3000.001, 3001.0):
        for ang in (0.0, 0.7, 2.0, 4.0):
            extra.append((r * math.cos(ang), 0.0, -7500.0 + r * math.sin(ang)))
    pts = sample_points(P["fault_box"], stations, extra)

    def field_cmp(name, mfn, mref, ssname, sref, xs=pts):
        worst, wp, wv = 0.0, None, None
        nbad = 0
        for x in xs:
            mv, sv = mfn(x), fault(ssname, x)
            if not close(mv, sv):
                nbad += 1
                d = abs(mv - sv) / max(abs(mv), abs(sv), 1e-300)
                if d > worst:
                    worst, wp, wv = d, x, (mv, sv)
        if nbad:
            desc = f"{nbad}/{len(xs)} pts differ; worst rel {worst:.3e} at {wp}: {wv[0]!r} vs {wv[1]!r}"
        else:
            desc = f"{len(xs)} pts agree"
        R.add(name, desc, mref, ssname, sref, nbad == 0)

    # ---- prestress tensor
    comp = {"s_xx": (0, 0), "s_yy": (1, 1), "s_zz": (2, 2), "s_xy": (0, 1), "s_yz": (1, 2), "s_xz": (0, 2)}
    sref_st = deck.ref("stress", "sigma_xy_pa") + (", [[stress_region]] x3" if "stress_region" in deck.t else "")
    for c, (i, j) in comp.items():
        field_cmp(f"prestress {c} [Pa]", lambda x, i=i, j=j: float(mfem_prestress(deck, x)[i, j]),
                  sref_st + " + driver:1567-1645", c, fy)

    # ---- nucleation
    nuc_names = {"s_xx": "nuc_xx", "s_yy": "nuc_yy", "s_zz": "nuc_zz", "s_xy": "nuc_xy",
                 "s_yz": "nuc_yz", "s_xz": "nuc_xz"}
    if "nucleation" in deck.t:
        nk = deck.t["nucleation"]["kind"]
        for c, (i, j) in comp.items():
            field_cmp(f"nucleation {nuc_names[c]} [Pa]",
                      lambda x, i=i, j=j: float(mfem_nucleation(deck, x)[i, j]),
                      deck.ref(f"nucleation.{nk}", "delta_tau_pa") + " + driver:1656-1684", nuc_names[c], fy)
        T = seconds(deck.get(f"nucleation.{nk}", "T_nuc_s"))
        t0 = par_num(par["dynamicrupture"]["t_0"][0])
        R.scalar("nucleation time T_nuc / t_0 [s]", T, deck.ref(f"nucleation.{nk}", "T_nuc_s"),
                 t0, pr("dynamicrupture", "t_0"))
        ts = [0.0, 1e-6, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999, 1.0, 1.5, 12.0]
        ok = all(close(mfem_smoothstep(t, T), seissol_smoothstep(t, t0)) for t in ts)
        R.add("G(t) at 11 times", "SmoothStep", "spatial_nucleation.cpp:20-29", "smoothStep",
              "GaussianNucleationFunction.h:22-31", ok)
    else:
        for c in comp:
            field_cmp(f"nucleation {nuc_names[c]} [Pa]", lambda x: 0.0, "no [nucleation] table",
                      nuc_names[c], fy)

    # ---- total SCEC tractions at t = 0, 0.5, 2 s, mapped through the to_scec.py rules
    pn, ps, pd = mfem_frame(deck)
    n_ss, s_ss, d_ss = seissol_frame(par)
    t0 = par_num(par["dynamicrupture"].get("t_0", ("0", 0))[0])
    Tm = seconds(deck.get(f"nucleation.{deck.t['nucleation']['kind']}", "T_nuc_s")) if "nucleation" in deck.t else 0.0
    bad = 0
    for x in pts:
        for t in (0.0, 0.5, 2.0):
            Sm = mfem_prestress(deck, x)
            dn = mfem_nucleation(deck, x)
            if dn is not None:
                Sm = Sm + dn * mfem_smoothstep(t, Tm)
            Tv = Sm @ pn
            m_tau1, m_tau2, m_sn = -float(Tv @ pd), -float(Tv @ ps), -float(Tv @ pn)
            Ss = np.zeros((3, 3))
            for c, (i, j) in comp.items():
                v = fault(c, x) + fault(nuc_names[c], x) * seissol_smoothstep(t, t0)
                Ss[i, j] = Ss[j, i] = v
            ts_, td_, pn_ = float(n_ss @ Ss @ s_ss), float(n_ss @ Ss @ d_ss), float(n_ss @ Ss @ n_ss)
            # to_scec.py: tau2 = Ts0, tau1 = Td0, sigma_n = -Pn0
            if not (close(m_tau2, ts_) and close(m_tau1, td_) and close(m_sn, -pn_)):
                bad += 1
    R.add("SCEC tau1/tau2/sigma_n (t=0,0.5,2)", "-T.dhat, -T.shat, -T.nhat",
          "driver:131-137, station_conventions", "Td0, Ts0, -Pn0 (n,strike,dip from par)",
          "ReceiverBasedOutput.cpp:308-325, Misc.cpp:19-34", bad == 0)
    vec = lambda v: "(" + ",".join(f"{float(c):g}" for c in v) + ")"
    R.add("frame strike / dip", f"s={vec(ps)} d={vec(pd)}", "driver:1647-1654",
          f"s={vec(s_ss)} d={vec(d_ss)}", "Misc.cpp:19-34",
          np.allclose(ps, s_ss) and np.allclose(pd, d_ss))

    # ---- friction
    if law == "slip_weakening":
        sw = deck.t["friction"]["slip_weakening"]
        for mk, sk in (("mu_s_default", "mu_s"), ("mu_d_default", "mu_d"), ("d_c_default", "d_c")):
            field_cmp(f"{sk}", lambda x, v=float(sw[mk]): v, deck.ref("friction.slip_weakening", mk), sk, fy)
        # SeisSol strength = -cohesion - mu sn (LinearSlipWeakening.h:158-160): cohesion_SS = -c_MFEM
        field_cmp("cohesion (SeisSol sign: -c)", lambda x, v=float(sw["cohesion_default"]): -v,
                  deck.ref("friction.slip_weakening", "cohesion_default"), "cohesion", fy)
        R.add("barrier / spatial rules", "none", deck.path.split("/")[-1], "none", fy,
              not sw.get("spatial") and "forced_rupture_time" not in fault.provides)
    else:
        rs = deck.t["friction"]["rate_state"]
        dr = par["dynamicrupture"]
        R.scalar("f0", rs["f_0_default"], deck.ref("friction.rate_state", "f_0_default"),
                 par_num(dr["rs_f0"][0]), pr("dynamicrupture", "rs_f0"))
        R.scalar("b", rs["b_default"], deck.ref("friction.rate_state", "b_default"),
                 par_num(dr["rs_b"][0]), pr("dynamicrupture", "rs_b"))
        R.scalar("V0 [m/s]", rs["V_0_default"], deck.ref("friction.rate_state", "V_0_default"),
                 par_num(dr["rs_sr0"][0]), pr("dynamicrupture", "rs_sr0"))
        vini = math.hypot(par_num(dr["rs_inisliprate1"][0]), par_num(dr["rs_inisliprate2"][0]))
        R.scalar("V_init [m/s] (magnitude)", rs["V_init_default"], deck.ref("friction.rate_state", "V_init_default"),
                 vini, pr("dynamicrupture", "rs_inisliprate1"))
        for k in ("rs_f0", "rs_b", "rs_muw"):
            if k in fault.provides:
                R.add(f"{k} in yaml", "-", "-", "present", fy, False)
        field_cmp("a", lambda x: mfem_rs_fields(deck, x)["a"],
                  deck.ref("friction.rate_state.spatial", "a_inner") + " boxcar_taper", "rs_a", fy)
        field_cmp("L = Dc [m]", lambda x: mfem_rs_fields(deck, x)["Dc"],
                  deck.ref("friction.rate_state", "Dc_default"), "rs_sl0", fy)
        if evo == "slip_law_strong_rate_weakening":
            R.scalar("f_w = muW", rs["f_w_default"], deck.ref("friction.rate_state", "f_w_default"),
                     par_num(dr["rs_muw"][0]), pr("dynamicrupture", "rs_muw"))
            field_cmp("V_w [m/s]", lambda x: mfem_rs_fields(deck, x)["V_w"],
                      deck.ref("friction.rate_state.spatial", "V_w_inner") + " boxcar_taper", "rs_srW", fy)
        # initial state from the t = 0 traction
        f0, bb, V0 = rs["f_0_default"], rs["b_default"], rs["V_0_default"]
        nbad, worst = 0, 0.0
        for x in pts:
            psi_m = mfem_psi0(deck, x)
            Ss = np.zeros((3, 3))
            for c, (i, j) in comp.items():
                Ss[i, j] = Ss[j, i] = fault(c, x)
            tr = Ss @ n_ss
            t1, t2, pnn = float(tr @ s_ss), float(tr @ d_ss), float(tr @ n_ss)
            aa, L = fault("rs_a", x), fault("rs_sl0", x)
            tmp = abs(math.hypot(t1, t2) / (aa * pnn))
            if fl == 103:
                psi_s = aa * logsinh(2.0 * V0 / vini, tmp)
            else:
                theta = L / V0 * math.exp((aa * logsinh(2.0, tmp) - f0 - aa * math.log(vini / V0)) / bb)
                psi_s = f0 + bb * math.log(V0 * theta / L)
            if not close(psi_m, psi_s):
                nbad += 1
                worst = max(worst, abs(psi_m - psi_s) / abs(psi_m))
        R.add("initial state psi_0", f"{len(pts) - nbad}/{len(pts)} agree" + (f", worst {worst:.2e}" if nbad else ""),
              "dieterich_ruina.hpp:507-529", "FL=103 psi_0" if fl == 103 else "FL=3 theta_0 -> psi",
              "RateAndStateInitializer.cpp:154-180" if fl == 103 else "RateAndStateInitializer.cpp:104-132",
              nbad == 0)

    # ---- not compared
    num = deck.t.get("numerics", {})
    R.info("Kelvin-Voigt damping", f"kv_eta={num.get('kv_eta')} kv_scale={num.get('kv_scale')}",
           deck.ref("numerics", "kv_eta"), "none (no KV in SeisSol a8e7f04f)", "-")
    if "prakash_clifton" in deck.t.get("friction", {}):
        pc = deck.t["friction"]["prakash_clifton"]
        R.info("Prakash-Clifton", f"v*={pc['v_star']} L={pc['length']} (on sigma_n)",
               deck.ref("friction.prakash_clifton", "v_star"),
               "none (only FL=6, on the strength: LinearSlipWeakening.h:280-295)", "-")
    R.info("CFL", f"cfl={num.get('cfl')} x C(p) (P3)", deck.ref("numerics", "cfl"),
           f"CFL={par['discretization']['cfl'][0]} ClusteredLTS={par['discretization']['clusteredlts'][0]}",
           pr("discretization", "cfl") + " (different definitions; speed, not physics)")
    R.show()
    print(f"\n{'PASS' if R.bad == 0 else 'FAIL'}: {len(R.rows) - R.bad}/{len(R.rows)} compared items match"
          f" (rel tol {REL:g}); {len(R.open)} open items")
    return 0 if R.bad == 0 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, KeyError, ValueError) as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(2)
