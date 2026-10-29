#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/rotate_case.py — 2D in-plane design rev 3, phase A: the
# BP5-mini rotation run (tests/2d_inplane/rotation_run.sh).  A rigid
# rotation R of the whole problem, the mesh and every vector of the deck,
# must leave every fault-local and axis-projected output unchanged.
#
#   rotate_case.py mesh IN.msh OUT.msh       the nodes of a Gmsh 2.2 ASCII mesh rotated by R
#   rotate_case.py deck IN.toml OUT.toml MESH
#                                            the deck with every vector rotated by R, [mesh].path
#                                            = MESH, and [cycle.loading] kind = "plate" with BP5's
#                                            plane, normal and slip direction rotated
#   rotate_case.py plate IN.toml OUT.toml    the unrotated reference: the deck with [cycle.loading]
#                                            kind = "plate" and BP5's vectors (equal to the default
#                                            "bp5" loading bit for bit, test_plate_loading L1)
#
# 2-D decks and meshes (the tensor prestress rotation run,
# tests/2d_inplane/tensor_rotation_run.sh; Chunhui Zhao, 2026-10-29: "test on
# the rotated planar fault line (modified from TPV205-2D) in a certain
# angle"):
#   rotate_case.py mesh2d IN.msh OUT.msh DEG  the nodes turned by DEG degrees about z
#   rotate_case.py deck2d IN.toml OUT.toml MESH DEG
#                                            every 2-vector of the deck turned by DEG
#                                            about z, [mesh].path = MESH, and [stress]'s
#                                            tensor S -> R S R^T (kind = "constant_tensor";
#                                            sigma_xz = sigma_yz = 0); the deck's
#                                            [cycle.loading] is turned like the rest
#
# R: 40 degrees about (1, 2, 3), as frame_box.hpp's Rot::Axis(1, 2, 3, 40).
# Rotated (points and directions alike, R is about the origin, and the
# deck's mesh_scale commutes with it): [fault_geometry] ref_normal, up;
# [cycle.fields.*] center, axis1, axis2; [cycle.output] plane_origin,
# plane_axis1, plane_axis2; [[cycle.output.station]] x; [cycle.benchmark]
# volume_axes (each row); [[cycle.benchmark.volume_station]] x;
# [cycle.loading] plane_point, plane_normal, slip_dir.  Not rotated: the
# 2-vectors of [cycle.benchmark] (region_x2, region_x3, profile from/to),
# which are coordinates in the output plane's own axes.  Any other array of
# three numbers is refused, so a new deck vector cannot be missed silently.
# deck2d: the same keys as 2-vectors; up must be [0, 0, +-1] (z, which the
# rotation keeps); any other array of two or three numbers is refused
# ([cycle.snapshots.volume] box_min and box_max included: an axis-aligned box
# does not turn).
#
# Scope (Chunhui Zhao, 2026-10-25).
import math
import re
import sys

try:
    import tomllib
except ImportError:      # Python < 3.11
    tomllib = None

AXIS, DEG = (1.0, 2.0, 3.0), 40.0
VECTOR_KEYS = {"ref_normal", "up", "center", "axis1", "axis2", "plane_origin", "plane_axis1",
               "plane_axis2", "x", "plane_point", "plane_normal", "slip_dir"}
MATRIX_KEYS = {"volume_axes"}
BP5_LOADING = {"plane_point": [0.0, 0.0, 0.0], "plane_normal": [0.0, 1.0, 0.0],
               "slip_dir": [1.0, 0.0, 0.0]}


def rotation(axis, deg):
    l = math.sqrt(sum(a * a for a in axis))
    k = [a / l for a in axis]
    t = math.radians(deg)
    c, s = math.cos(t), math.sin(t)
    R = [[(c if i == j else 0.0) + (1 - c) * k[i] * k[j] for j in range(3)] for i in range(3)]
    R[0][1] -= s * k[2]; R[0][2] += s * k[1]
    R[1][0] += s * k[2]; R[1][2] -= s * k[0]
    R[2][0] -= s * k[1]; R[2][1] += s * k[0]
    return R


def apply(R, v):
    return [R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2] for i in range(3)]


def fmt(v):
    return "[" + ", ".join(repr(float(a)) for a in v) + "]"


def parse_value(text):
    if tomllib is None:
        return eval(text, {"__builtins__": {}}, {"true": True, "false": False})
    return tomllib.loads("v = " + text)["v"]


def is_vec3(v):
    return isinstance(v, list) and len(v) == 3 and all(isinstance(a, (int, float)) for a in v)


def rotate_deck(src, dst, R, mesh_path):
    out, section = [], ""
    line_re = re.compile(r"^(\s*)([A-Za-z0-9_]+)(\s*=\s*)(\[.*\])(\s*(#.*)?)$")
    for line in open(src).read().split("\n"):
        h = re.match(r"^\s*\[\[?([^\]]+)\]\]?", line)
        if h and "=" not in line.split("#")[0]:
            section = h.group(1).strip()
        m = line_re.match(line)
        if m:
            key, val = m.group(2), parse_value(m.group(4).split("#")[0].strip())
            if key in MATRIX_KEYS:
                assert isinstance(val, list) and all(is_vec3(r) for r in val), (section, key)
                line = m.group(1) + key + m.group(3) + "[" + ", ".join(fmt(apply(R, r)) for r in val) + "]" + m.group(5)
            elif key in VECTOR_KEYS:
                assert is_vec3(val), (section, key, val)
                line = m.group(1) + key + m.group(3) + fmt(apply(R, val)) + m.group(5)
            elif is_vec3(val):
                sys.exit(f"rotate_case.py: [{section}] {key} = {val} is an array of three numbers "
                         "that is not in the list of deck vectors: refusing")
        if mesh_path and section == "mesh" and re.match(r"^\s*path\s*=", line):
            line = f'path  = "{mesh_path}"'
        out.append(line)
    text = "\n".join(out)
    assert "[cycle.loading]" not in text, "the deck already has a [cycle.loading] table"
    L = {k: apply(R, v) for k, v in BP5_LOADING.items()}
    text = text.rstrip("\n") + "\n\n[cycle.loading]\nkind         = \"plate\"\n" + \
        "".join(f"{k:<12} = {fmt(L[k])}\n" for k in ("plane_point", "plane_normal", "slip_dir")) + \
        "band         = 1000.0\n"
    open(dst, "w").write(text)


PLANE_KEYS = {"region_x2", "region_x3", "from", "to"}   # output-plane coordinates
TENSOR_KEYS = ("sigma_xx_pa", "sigma_yy_pa", "sigma_zz_pa", "sigma_xy_pa", "sigma_yz_pa",
               "sigma_xz_pa")


def rotation_z(deg):
    return rotation((0.0, 0.0, 1.0), deg)


def is_vec(v, n):
    return isinstance(v, list) and len(v) == n and all(isinstance(a, (int, float)) for a in v)


def rotate_deck_2d(src, dst, deg, mesh_path):
    R = rotation_z(deg)
    text = open(src).read()
    deck = tomllib.loads(text)
    st = deck.get("stress", {})
    rot_tensor = {}
    if st.get("kind") == "constant_tensor":
        g = {k: float(st[k]) for k in TENSOR_KEYS}
        assert g["sigma_xz_pa"] == 0.0 and g["sigma_yz_pa"] == 0.0, "a 2-D deck with anti-plane stress"
        S = [[g["sigma_xx_pa"], g["sigma_xy_pa"], 0.0], [g["sigma_xy_pa"], g["sigma_yy_pa"], 0.0],
             [0.0, 0.0, g["sigma_zz_pa"]]]
        RS = [[sum(R[i][k] * S[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        Sp = [[sum(RS[i][k] * R[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
        rot_tensor = {"sigma_xx_pa": Sp[0][0], "sigma_yy_pa": Sp[1][1], "sigma_zz_pa": Sp[2][2],
                      "sigma_xy_pa": 0.5 * (Sp[0][1] + Sp[1][0]), "sigma_yz_pa": 0.0,
                      "sigma_xz_pa": 0.0}
    out, section, done = [], "", set()
    arr_re = re.compile(r"^(\s*)([A-Za-z0-9_]+)(\s*=\s*)(\[.*\])(\s*(#.*)?)$")
    num_re = re.compile(r"^(\s*)(sigma_[xyz]{2}_pa)(\s*=\s*)([^#\s]+)(\s*(#.*)?)$")
    rot2 = lambda v: [R[0][0] * v[0] + R[0][1] * v[1], R[1][0] * v[0] + R[1][1] * v[1]]
    for line in text.split("\n"):
        h = re.match(r"^\s*\[\[?([^\]]+)\]\]?", line)
        if h and "=" not in line.split("#")[0]:
            section = h.group(1).strip()
        m = arr_re.match(line)
        if m:
            key, val = m.group(2), parse_value(m.group(4).split("#")[0].strip())
            if key in MATRIX_KEYS:
                assert isinstance(val, list) and all(is_vec(r, 2) for r in val), (section, key)
                line = m.group(1) + key + m.group(3) + "[" + ", ".join(fmt(rot2(r)) for r in val) + "]" + m.group(5)
            elif key == "up":
                assert is_vec(val, 3) and val[0] == 0 and val[1] == 0 and abs(val[2]) == 1, \
                    f"[{section}] up = {val}: deck2d needs up = [0, 0, +-1]"
            elif key in VECTOR_KEYS:
                assert is_vec(val, 2), (section, key, val)
                line = m.group(1) + key + m.group(3) + fmt(rot2(val)) + m.group(5)
            elif key in PLANE_KEYS and section.startswith("cycle.benchmark"):
                assert is_vec(val, 2), (section, key, val)
            elif is_vec(val, 2) or is_vec(val, 3):
                sys.exit(f"rotate_case.py: [{section}] {key} = {val} is an array of numbers "
                         "that is not in the list of deck vectors: refusing")
        n = num_re.match(line)
        if n and section == "stress" and rot_tensor:
            key = n.group(2)
            line = n.group(1) + key + n.group(3) + repr(float(rot_tensor[key])) + n.group(5)
            done.add(key)
        if mesh_path and section == "mesh" and re.match(r"^\s*path\s*=", line):
            line = f'path  = "{mesh_path}"'
        out.append(line)
    assert not rot_tensor or done == set(TENSOR_KEYS), done
    open(dst, "w").write("\n".join(out))


def rotate_mesh(src, dst, R):
    lines = open(src).read().split("\n")
    assert lines[0] == "$MeshFormat" and lines[1].split()[0] == "2.2", "Gmsh 2.2 ASCII only"
    i = lines.index("$Nodes")
    n = int(lines[i + 1])
    for j in range(i + 2, i + 2 + n):
        f = lines[j].split()
        x = apply(R, [float(f[1]), float(f[2]), float(f[3])])
        lines[j] = f[0] + " " + " ".join(repr(a) for a in x)
    assert lines[i + 2 + n] == "$EndNodes"
    open(dst, "w").write("\n".join(lines))


def main(argv):
    mode = argv[1] if len(argv) > 1 else ""
    R = rotation(AXIS, DEG)
    I = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    if mode == "mesh" and len(argv) == 4:
        rotate_mesh(argv[2], argv[3], R)
    elif mode == "deck" and len(argv) == 5:
        rotate_deck(argv[2], argv[3], R, argv[4])
    elif mode == "plate" and len(argv) == 4:
        rotate_deck(argv[2], argv[3], I, None)
    elif mode == "mesh2d" and len(argv) == 5:
        rotate_mesh(argv[2], argv[3], rotation_z(float(argv[4])))
    elif mode == "deck2d" and len(argv) == 6:
        rotate_deck_2d(argv[2], argv[3], float(argv[5]), argv[4])
    else:
        sys.exit(__doc__ or "usage: rotate_case.py mesh IN OUT | deck IN OUT MESH | plate IN OUT"
                 " | mesh2d IN OUT DEG | deck2d IN OUT MESH DEG")


if __name__ == "__main__":
    main(sys.argv)
