#!/usr/bin/env python3
"""apply_alt_patches.py -- de-PREFERRED the copied toolchain.

The ALT code/ tree is a byte-for-byte copy of the PREFERRED ShakeOut-D build, so
every PREFERRED-specific constant in it is now silently wrong for ALT. Each edit
below was verified against the file before being written, and each is a value
that would NOT raise -- it would produce a plausible wrong answer.

  1. fill_to_msh.py hardcodes PREFERRED's fault facet count (2,761,488) and
     fault area (7.96255e9 m2). ALT has 2,564,480 / ~6.5376e9. Left alone the
     pinhole count and the hull-count assert are both nonsense.
  2. gate_census2.py / relax_rv.py / make_view_xdmf.py / leb_gate_close.py carry
     BOX = (72000, 786000, 3524000, 3996000) -- the OLD, DELETED axis-aligned
     ShakeOut box. It is 13.4 km short on the south and 45.7 km short on the
     north of ShakeOut-D. material.py clips with np.searchsorted, so those wedges
     silently sample the EDGE ROW of the MUSCAL subset: no NaN, no warning, and
     the coverage counter still reads 100 %. This is wrong for PREFERRED too.
  3. gate_census2.py computes its diagnostics under `if k == "both"`, i.e. on
     min(deck, MUSCAL). The spec says score on native MUSCAL, never the deck
     cube -- the deck's 250 m z-binning is ~6x more demanding.

ZTOP is deliberately NOT patched here: it is a design decision (see OQ-1), not a
defect.
"""
import re
import sys
from pathlib import Path

# ShakeOut-D bounding box in UTM 11N, from domain_corners_utm.npy:
#   x 132,679.3 .. 723,873.2   y 3,490,617.1 .. 4,054,679.7
NEW_BOX = "(132000.0, 724000.0, 3490000.0, 4055000.0)"
NEW_BOX_LIST = "[132000.0, 724000.0, 3490000.0, 4055000.0]"
OLD_BOX_T = "(72000.0, 786000.0, 3524000.0, 3996000.0)"
OLD_BOX_L = "[72000.0, 786000.0, 3524000.0, 3996000.0]"


def patch(path, subs, must=True):
    p = Path(path)
    s0 = p.read_text()
    s = s0
    for old, new in subs:
        if old not in s:
            if must:
                print(f"  !! NOT FOUND in {p.name}: {old[:70]}")
            continue
        s = s.replace(old, new)
    if s != s0:
        p.write_text(s)
        print(f"  patched {p.name}")
    else:
        print(f"  UNCHANGED {p.name}")


def main():
    code = Path(sys.argv[1])
    print("[1] fill_to_msh.py -- PREFERRED fault count / area -> derived")
    patch(code / "fill_to_msh.py", [
        ("miss = 2761488 - nuniq",
         "NF = int((MARK == 7).sum())   # ALT: derive, never hardcode a lineage's count\n"
         "miss = NF - nuniq"),
        ("== len(plcT) - 2761488", "== len(plcT) - NF"),
        ("100*ar.sum()/7.96255e9", "100*ar.sum()/FAULT_AREA"),
    ])
    # define FAULT_AREA next to NF
    p = code / "fill_to_msh.py"
    s = p.read_text()
    if "FAULT_AREA" in s and "FAULT_AREA =" not in s:
        s = s.replace(
            "NF = int((MARK == 7).sum())",
            "NF = int((MARK == 7).sum())\n"
            "_fp = plcP[plcT[MARK == 7]]\n"
            "FAULT_AREA = float(0.5*np.linalg.norm(np.cross(_fp[:,1]-_fp[:,0], _fp[:,2]-_fp[:,0]), axis=1).sum())")
        p.write_text(s)
        print("  added FAULT_AREA derivation")

    print("[2] stale ShakeOut BOX -> ShakeOut-D")
    for f in ("gate_census2.py", "relax_rv.py", "make_view_xdmf.py"):
        patch(code / f, [(OLD_BOX_T, NEW_BOX)])
    patch(code / "leb_gate_close.py", [(OLD_BOX_L, NEW_BOX_LIST)])

    print("[3] gate_census2.py -- diagnostics on native MUSCAL, not min(deck,native)")
    patch(code / "gate_census2.py", [('if k == "both":', 'if k == "native":')])

    print("\n[verify] no PREFERRED constants left:")
    bad = 0
    for f in sorted(code.glob("*.py")):
        t = f.read_text()
        for pat in ("2761488", "7.96255e9", "72000.0, 786000.0"):
            if pat in t:
                print(f"  !! {f.name} still contains {pat}")
                bad += 1
    print("  clean" if not bad else f"  {bad} remaining")


if __name__ == "__main__":
    main()
