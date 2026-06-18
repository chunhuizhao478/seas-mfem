#!/usr/bin/env python3
"""Phase 2 verification: check the pumgen-produced PUML for the deep multi-strand
mesh (PLAN_safv4_deep_mesh_port_2026-06-15.md).

Run AFTER `bash build_and_run_pumgen_linux.sh .../safs_mesh_deep.msh` produces
`safs_mesh_deep.puml.h5`.  Runnable on Linux (in the pumgen-build conda env, which
has h5py+numpy) or back on the Mac (pythonenv).

Checks, against the Phase-0/1 verified `.msh` triangle counts:
  - mesh size: connect (tets) == 1,180,159 ; geometry (nodes) == 279,448
  - boundary codes are a subset of {0 interior, 1 free surface, 3 DR, 5 absorbing}
  - free surface (BC 1) face count == 175,033   (exact; true boundary)
  - absorbing   (BC 5) face count == 122,593    (exact; true boundary)
  - dynamic rupture (BC 3) == 132,333 (1x) OR 264,666 (2x)  -- DR is an INTERNAL
    surface; some pumgen versions tag both adjacent tet faces, so 1x or 2x the
    132,333 fault triangles are both acceptable.  Reported, not hard-failed beyond
    that set.

Usage: verify_puml_deep.py [safs_mesh_deep.puml.h5]
Exit 0 on PASS, non-zero on FAIL.
"""
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("FAIL: h5py not available -- run in the pumgen-build conda env "
             "(Linux) or pythonenv (Mac).")

EXPECT_TETS = 1180159
EXPECT_NODES = 279448
EXPECT_FS = 175033        # BC 1 free surface  (= former tag 201)
EXPECT_ABS = 122593       # BC 5 absorbing     (= former tags 202+203)
EXPECT_DR_1X = 132333     # BC 3 dynamic rupture (= former tags 101+102+103)
EXPECT_DR_2X = 2 * EXPECT_DR_1X
ALLOWED_CODES = {0, 1, 3, 5}

# the verified PUML belongs in the parent (submission) folder, next to parameters.par
DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "safs_mesh_deep.puml.h5")
BC_NAME = {0: "interior", 1: "free surface", 3: "dynamic rupture", 5: "absorbing"}


def decode_boundary(b):
    """Return {bc_code: face_count}.

    Handles both PUML boundary layouts:
      - packed: shape (nElem,), 4 faces x 8 bits each
      - unpacked: shape (nElem, 4), one BC code per face
    """
    faces = {}
    if b.ndim == 2 and b.shape[1] == 4:
        codes, counts = np.unique(b.astype(np.int64), return_counts=True)
        for c, n in zip(codes, counts):
            faces[int(c)] = faces.get(int(c), 0) + int(n)
        return faces, "unpacked (nElem,4)"
    bb = b.astype(np.int64).ravel()
    for k in range(4):
        codes, counts = np.unique((bb >> (8 * k)) & 0xFF, return_counts=True)
        for c, n in zip(codes, counts):
            faces[int(c)] = faces.get(int(c), 0) + int(n)
    return faces, "packed 8-bit x4"


def main(argv):
    path = os.path.abspath(argv[1] if len(argv) > 1 else DEFAULT)
    print(f"puml: {path}")
    if not os.path.isfile(path):
        print(f"FAIL: not found: {path}")
        return 2

    fails = []
    with h5py.File(path, "r") as f:
        if "connect" not in f or "geometry" not in f or "boundary" not in f:
            print("FAIL: missing connect/geometry/boundary dataset; keys =",
                  list(f.keys()))
            return 2
        n_tets = f["connect"].shape[0]
        n_nodes = f["geometry"].shape[0]
        b = f["boundary"][:]

    print(f"connect (tets):  {n_tets}   (expect {EXPECT_TETS})")
    print(f"geometry (nodes):{n_nodes}   (expect {EXPECT_NODES})")
    if n_tets != EXPECT_TETS:
        fails.append(f"tet count {n_tets} != {EXPECT_TETS}")
    if n_nodes != EXPECT_NODES:
        fails.append(f"node count {n_nodes} != {EXPECT_NODES}")

    faces, layout = decode_boundary(b)
    print(f"boundary layout: {layout}")
    print("boundary face counts by BC code:")
    for c in sorted(faces):
        print(f"   BC {c} ({BC_NAME.get(c, '?')}): {faces[c]}")

    nonzero = {c: n for c, n in faces.items() if c != 0}
    if not set(nonzero).issubset(ALLOWED_CODES):
        fails.append(f"boundary codes {set(nonzero)} not subset of {ALLOWED_CODES}")

    fs = faces.get(1, 0)
    ab = faces.get(5, 0)
    dr = faces.get(3, 0)
    if fs != EXPECT_FS:
        fails.append(f"free-surface (BC1) {fs} != {EXPECT_FS}")
    if ab != EXPECT_ABS:
        fails.append(f"absorbing (BC5) {ab} != {EXPECT_ABS}")
    if dr == EXPECT_DR_1X:
        print(f"DR (BC3) = {dr} = 1x fault triangles (single-sided tagging)")
    elif dr == EXPECT_DR_2X:
        print(f"DR (BC3) = {dr} = 2x fault triangles (both adjacent tets tagged)")
    else:
        fails.append(f"DR (BC3) {dr} != {EXPECT_DR_1X} (1x) and != "
                     f"{EXPECT_DR_2X} (2x)")

    print("\n" + "=" * 60)
    if fails:
        print("PUML VERIFY: FAIL")
        for x in fails:
            print("  - " + x)
        return 1
    print("PUML VERIFY: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
