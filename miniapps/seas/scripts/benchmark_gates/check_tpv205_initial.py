#!/usr/bin/env python3
"""Gates D-01 and D-02 (benchmark plan rev 1.3) on <prefix>_initial_r*.csv of a TPV205 Nitsche run.

D-01: at every fault DOF, tau0 = -T.shat is the value of the face the DOF belongs to
      (70, 81.6, 78.0 or 62.0 MPa) and sigma_n0 = 120 MPa, both to round-off.  A DOF strictly
      inside a patch must carry the patch value, one strictly outside every patch 70 MPa; a DOF on
      a patch edge carries either (one DG copy per adjacent face).
D-02: mu_s, mu_d, Dc, cohesion equal the deck values at every DOF; zero DOFs carry a barrier mu_s.

Usage: check_tpv205_initial.py <output_dir> <prefix>
"""
import glob, sys
out, prefix = sys.argv[1], sys.argv[2]
rows = []
for f in sorted(glob.glob(f"{out}/{prefix}_initial_r*.csv")):
    hdr = None
    for line in open(f):
        if line.startswith('#'):
            if hdr is None: hdr = line[1:].split()
            continue
        rows.append([float(v) for v in line.split()])
assert hdr == ['x','y','z','a','Dc','psi0','tau0_strike(-T.shat)','sigma_n0','V0_strike','mu_s','mu_d','cohesion'], hdr
col = {k: i for i, k in enumerate(hdr)}
patches = [(-1500, 1500, 81.6e6), (-9000, -6000, 78.0e6), (6000, 9000, 62.0e6)]
zlo, zhi = -9000.0, -6000.0
tol_pos = 1e-3          # m: "on an edge" band
tol_rel = 1e-13
nin = nout = nedge = 0; bad = []
maxdev_tau = maxdev_sig = 0.0
for r in rows:
    x, z = r[col['x']], r[col['z']]
    tau, sig = r[col['tau0_strike(-T.shat)']], r[col['sigma_n0']]
    maxdev_sig = max(maxdev_sig, abs(sig - 120e6) / 120e6)
    status, val = 'out', None
    for (a, b, v) in patches:
        strictly_in = a + tol_pos < x < b - tol_pos and zlo + tol_pos < z < zhi - tol_pos
        near = a - tol_pos <= x <= b + tol_pos and zlo - tol_pos <= z <= zhi + tol_pos
        if strictly_in: status, val = 'in', v
        elif near and status != 'in': status, val = 'edge', v
    if status == 'in':
        nin += 1; ok = abs(tau - val) <= tol_rel * val; dev = abs(tau - val) / val
    elif status == 'out':
        nout += 1; ok = abs(tau - 70e6) <= tol_rel * 70e6; dev = abs(tau - 70e6) / 70e6
    else:
        nedge += 1
        dev = min(abs(tau - val) / val, abs(tau - 70e6) / 70e6); ok = dev <= tol_rel
    maxdev_tau = max(maxdev_tau, dev)
    if not ok: bad.append((x, z, tau, status, val))
print(f"D-01: {len(rows)} DOFs: {nin} strictly inside a patch, {nout} strictly outside, {nedge} on a patch edge; "
      f"max relative deviation of tau0 {maxdev_tau:.3g}, of sigma_n0 {maxdev_sig:.3g}")
for b in bad[:10]: print("   BAD", b)
d01 = not bad and maxdev_sig <= tol_rel and nin > 0 and nedge > 0
print("  PASSED: D-01" if d01 else "  FAILED: D-01")
deck = {'mu_s': 0.677, 'mu_d': 0.525, 'Dc': 0.40, 'cohesion': 0.0}
d02 = True
for k, v in deck.items():
    vals = {r[col[k]] for r in rows}
    ok = vals == {v}
    d02 &= ok
    print(f"D-02: {k} values {sorted(vals)[:4]} (deck {v}) {'ok' if ok else 'MISMATCH'}")
nbar = sum(1 for r in rows if r[col['mu_s']] >= 1e6)
print(f"D-02: DOFs with a barrier mu_s: {nbar}")
d02 &= nbar == 0
print("  PASSED: D-02" if d02 else "  FAILED: D-02")
sys.exit(0 if d01 and d02 else 1)
