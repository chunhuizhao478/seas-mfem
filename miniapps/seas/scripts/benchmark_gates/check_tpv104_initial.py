#!/usr/bin/env python3
"""Gate D-06 (benchmark plan rev 1.3) and the per-node part of S-02 on
<prefix>_initial_r*.csv of a TPV104 Nitsche run.

D-06: a and V_w at every fault DOF equal the SCEC TPV104 eqs. 8-10 (smooth boxcar, a_in 0.01,
      a_out 0.02, V_w_in 0.1, V_w_out 1.0, Ls 15 km, W/2 7.5 km about 7.5 km depth, w 3 km)
      to 1e-15 x max(1, |value|).
S-02: psi0 = a ln((2 V0 / V_ini) sinh(tau0 / (sigma_n a))) at every DOF (eq. 11), 1e-14 relative,
      with tau0 and sigma_n0 read from the same row.

Usage: check_tpv104_initial.py <output_dir> <prefix>
"""
import glob, math, sys
out, prefix = sys.argv[1], sys.argv[2]
rows, hdr = [], None
for f in sorted(glob.glob(f"{out}/{prefix}_initial_r*.csv")):
    h = None
    for line in open(f):
        if line.startswith('#'):
            if h is None: h = line[1:].split()
            continue
        rows.append([float(v) for v in line.split()])
    hdr = h
assert hdr == ['x','y','z','a','Dc','psi0','tau0_strike(-T.shat)','sigma_n0','V0_strike','V_w'], hdr
c = {k: i for i, k in enumerate(hdr)}
def boxcar(x, W, w):
    ax = abs(x)
    if ax <= W: return 1.0
    if ax >= W + w: return 0.0
    return 0.5 * (1.0 + math.tanh(w / (ax - W - w) + w / (ax - W)))
da = dv = dp = 0.0
for r in rows:
    x, z = r[c['x']], r[c['z']]
    B = boxcar(x, 15e3, 3e3) * boxcar(-z - 7.5e3, 7.5e3, 3e3)
    a = 0.02 + (0.01 - 0.02) * B
    vw = 1.0 + (0.1 - 1.0) * B
    da = max(da, abs(r[c['a']] - a) / max(1.0, abs(a)))
    dv = max(dv, abs(r[c['V_w']] - vw) / max(1.0, abs(vw)))
    aa, tau, sig = r[c['a']], r[c['tau0_strike(-T.shat)']], r[c['sigma_n0']]
    psi = aa * math.log((2 * 1e-6 / 1e-16) * math.sinh(tau / (sig * aa)))
    dp = max(dp, abs(r[c['psi0']] - psi) / abs(psi))
amin = min(r[c['a']] for r in rows); amax = max(r[c['a']] for r in rows)
print(f"D-06: {len(rows)} DOFs, a in [{amin}, {amax}]; max |a - eq.| {da:.3g}, max |V_w - eq.| {dv:.3g}")
print(f"S-02 (per node): max relative |psi0 - eq. 11| {dp:.3g}")
ok6 = da <= 1e-15 and dv <= 1e-15
ok2 = dp <= 1e-14
print("  PASSED: D-06" if ok6 else "  FAILED: D-06")
print("  PASSED: S-02 per node" if ok2 else "  FAILED: S-02 per node")
sys.exit(0 if ok6 and ok2 else 1)
