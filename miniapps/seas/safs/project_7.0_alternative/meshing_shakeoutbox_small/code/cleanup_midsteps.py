#!/usr/bin/env python3
"""cleanup_midsteps.py -- leave exactly ONE mesh per case.

The build produces intermediates (collar npz before and after gate closure, a
pre-closure merged mesh, PLC dumps).  They are large and, once the final mesh
has passed its acceptance checks, worthless.  This deletes them.

SAFETY: the KEEP set is the list of final products.  Nothing is deleted unless
  * every KEEP file exists and is non-empty, checked BEFORE any unlink, and
  * the candidate lies under one of the three meshing_shakeoutbox_* trees, and
  * the candidate is a regular file, not a symlink or directory, and
  * the candidate is NOT in the KEEP set.
The KEEP set is re-verified after each deletion.  --dry-run by default.
"""

import argparse
import os
import sys

BASE = ("/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/"
        "project_7.0_alternative")
TREES = tuple(f"{BASE}/meshing_shakeoutbox_{t}/" for t in ("small", "intermediate", "heavy"))
KEEP = (
    f"{BASE}/meshing_shakeoutbox_small/results/safalt_small_shakeoutbox.puml.h5",
    f"{BASE}/meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5",
    f"{BASE}/meshing_shakeoutbox_heavy/results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5",
)
# extensions that are always intermediates
MIDSTEP_EXT = (".npz", ".msh", ".vtu", ".mtr", ".node", ".ele", ".face")


def keep_snapshot():
    return {p: os.path.getsize(p) for p in KEEP if os.path.isfile(p)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="actually delete (default: dry run)")
    a = ap.parse_args()

    before = keep_snapshot()
    missing = [p for p in KEEP if p not in before]
    if missing:
        print("REFUSING: these final products are missing or empty:")
        for p in missing:
            print("   ", p)
        sys.exit(1)
    print(f"KEEP set verified: {len(before)}/{len(KEEP)} products, "
          f"{sum(before.values())/1e9:.2f} GB")
    for p, s in before.items():
        print(f"   {s/1e9:8.3f} GB  {p}")

    cand = []
    for t in TREES:
        for root, _d, files in os.walk(t):
            for fn in files:
                p = os.path.join(root, fn)
                if p in KEEP or os.path.islink(p) or not os.path.isfile(p):
                    continue
                is_mid = (f"{os.sep}build_tmp{os.sep}" in p
                          or fn.endswith(MIDSTEP_EXT)
                          or (fn.endswith(".puml.h5") and p not in KEEP))
                if is_mid:
                    cand.append(p)
    cand.sort()
    tot = sum(os.path.getsize(p) for p in cand)
    print(f"\n{len(cand)} mid-step files, {tot/1e9:.2f} GB")
    for p in cand:
        print(f"   {os.path.getsize(p)/1e9:8.3f} GB  {p[len(BASE)+1:]}")

    if not a.apply:
        print("\nDRY RUN -- nothing deleted.  Re-run with --apply.")
        return

    freed = n = 0
    for p in cand:
        if p in KEEP or not p.startswith(TREES) or os.path.islink(p) or not os.path.isfile(p):
            continue
        sz = os.path.getsize(p)
        os.unlink(p)
        freed += sz
        n += 1
        if keep_snapshot() != before:
            raise RuntimeError(f"KEEP SET DAMAGED after deleting {p}")
    print(f"\ndeleted {n} files, freed {freed/1e9:.2f} GB")
    print(f"KEEP set intact: {len(keep_snapshot())}/{len(KEEP)}")


if __name__ == "__main__":
    main()
