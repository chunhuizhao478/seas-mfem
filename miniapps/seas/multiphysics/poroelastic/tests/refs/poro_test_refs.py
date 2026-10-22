#!/usr/bin/env python3
"""Reference values for the poroelastic unit and verification tests (plan rev 3, sections 6-7).

Every number a test compares against is computed here from the definitions of the derivation
document (rev 1.2): eq. (4.4) for the constants, (4.6) for the high-frequency speeds, (8.4) for
the relaxation time, (9.3)-(9.4) for the fitted coefficient, section 7.4 for the tpv5-poro start,
and the two patch-test fields PT-1 and PT-2 of plan section 6.3.  It replaces biot.py of plan
rev 2.2, which used K_s = 60 GPa for the TPV205 set while the plan's table used 100 GPa.

Usage: python3 poro_test_refs.py [--json refs.json] [--header poro_test_refs.hpp]
       python3 poro_test_refs.py --check DIR   (make check-poro-refs: regenerate refs.json and
                                                the header and compare them with DIR)
--check fails on any byte of the header that differs and on any refs.json value that differs,
except in the eigenvalue results modes_1km and dispersion_tpv5 (EIGEN_KEYS).  Those carry the
eigensolver's round-off, about machine epsilon times the largest eigenvalue of each solve
(3.66e5 1/s for tpv5-poro), so their last digits depend on the BLAS build; their differences
are printed and do not fail.  Every value the C++ tests read is in the header.
Chunhui Zhao's project, 2026-10-20.  Inputs: SeisSol tpv5-poro material.yaml (precomputed-seissol),
the TPV205 undrained-equivalent set of plan section 3.2, and the V-F3 contrast set (step P1.1):
the same undrained moduli and density as the TPV205-equivalent set, with porosity 0.2 and
permeability 1e-13 m^2 (100 times larger), so a different Skempton coefficient.
"""
import argparse
import filecmp
import json
import math
import os
import sys
import tempfile

import numpy as np
from scipy.optimize import brentq


def derive(mu, lam, Kf, Ks, kap, phi, rhof, rhos, T, eta):
    """Eq. (4.4): the derived constants, with SeisSol Helper.h's formulas."""
    Kd = lam + 2.0 / 3.0 * mu
    alpha = 1.0 - Kd / Ks
    M = Ks / (1.0 - phi - Kd / Ks + phi * Ks / Kf)          # Helper.h form
    M2 = 1.0 / (phi / Kf + (alpha - phi) / Ks)             # textbook form, same value
    lamu = lam + alpha ** 2 * M
    Ku = Kd + alpha ** 2 * M
    rho = (1.0 - phi) * rhos + phi * rhof
    m = T * rhof / phi
    b = eta / kap
    tau = (rho * m - rhof ** 2) / (b * rho)                 # eq. (8.4)
    Mr = np.array([[rho, rhof], [rhof, m]])
    Dr = np.array([[0.0, 0.0], [0.0, b]])
    ev_relax = np.linalg.eigvals(np.linalg.solve(Mr, Dr)).real.max()
    Kp = np.array([[lamu + 2 * mu, alpha * M], [alpha * M, M]])
    ev = np.sort(np.linalg.eigvals(np.linalg.solve(Mr, Kp)).real)   # eq. (4.6)
    d = dict(mu=mu, lam=lam, Kf=Kf, Ks=Ks, kappa=kap, phi=phi, rho_f=rhof, rho_s=rhos, T=T, eta=eta,
             Kd=Kd, alpha=alpha, M=M, M_textbook_rel_diff=abs(M - M2) / M, lam_u=lamu, Ku=Ku,
             B=alpha * M / Ku, rho=rho, m=m, b=b, tau=tau, tau_from_eig=1.0 / ev_relax,
             rho_m_minus_rhof2=rho * m - rhof ** 2,
             f_c=b / (2 * math.pi * m),                      # eta phi / (2 pi T rho_f kappa)
             c_fast=math.sqrt(ev[1]), c_slow=math.sqrt(ev[0]),
             c_s_hf=math.sqrt(mu / (rho - rhof ** 2 / m)), c_s_lf=math.sqrt(mu / rho),
             c_p_undrained=math.sqrt((lamu + 2 * mu) / rho), c_p_drained=math.sqrt((lam + 2 * mu) / rho),
             diffusivity=kap / eta * M * (lam + 2 * mu) / (lamu + 2 * mu))
    return d


def tpv205_equivalent(Ks=100e9, phi=0.1, kap=1e-15):
    """Plan section 3.2: undrained moduli and bulk density equal TPV205's.

    phi = 0.2, kap = 1e-13 gives the V-F3 contrast set (same undrained moduli and density)."""
    rho, cs, cp = 2670.0, 3464.0, 6000.0
    mu = rho * cs ** 2
    lamu_target = rho * (cp ** 2 - 2 * cs ** 2)
    rhof, Kf = 1000.0, 2.2e9
    Ku = lamu_target + 2.0 / 3.0 * mu

    def f(Kd):
        a = 1.0 - Kd / Ks
        M = 1.0 / (phi / Kf + (a - phi) / Ks)
        return Kd + a * a * M - Ku

    Kd = brentq(f, 1e9, Ks * 0.999, xtol=1e-6, rtol=1e-15)
    lam = Kd - 2.0 / 3.0 * mu
    rhos = (rho - phi * rhof) / (1 - phi)
    d = derive(mu, lam, Kf, Ks, kap, phi, rhof, rhos, 2.0, 1e-3)
    d['lam_u_target'] = lamu_target
    d['lam_u_rel_err'] = abs(d['lam_u'] - lamu_target) / lamu_target
    return d


def dispersion(d, freqs):
    """Complex wavenumbers of the fast P, slow P and S plane waves with drag, e^{i(kx - wt)}."""
    Mr = np.array([[d['rho'], d['rho_f']], [d['rho_f'], d['m']]])
    Kp = np.array([[d['lam_u'] + 2 * d['mu'], d['alpha'] * d['M']], [d['alpha'] * d['M'], d['M']]])
    out = []
    for f in freqs:
        w = 2 * math.pi * f
        A = w ** 2 * Mr + 1j * w * np.array([[0, 0], [0, d['b']]])   # K k^2 X = (w^2 M + i w D) X
        k2 = np.linalg.eigvals(np.linalg.solve(Kp, A))
        ks = sorted([np.sqrt(z) for z in k2], key=lambda k: k.real)   # fast P has the smaller Re k
        kS = np.sqrt(w ** 2 / d['mu'] * (d['rho'] - d['rho_f'] ** 2 * w / (d['m'] * w + 1j * d['b'])))
        row = dict(f=f)
        for name, k in (('fastP', ks[0]), ('slowP', ks[1]), ('S', kS)):
            row[name + '_speed'] = w / k.real
            row[name + '_Q'] = k.real / (2 * abs(k.imag)) if k.imag != 0 else float('inf')
        out.append(row)
    return out


def modes_real_k(d, wavelength):
    """Complex omega of the P and S modes for a real wavenumber (periodic-box plane waves, V-B1).

    P: det(K k^2 - w^2 M_rho - i w D) = 0, solved as a 4 x 4 linear eigenproblem in w;
    S: mu k^2 X_u = (w^2 M_rho + i w D) restricted to the shear rows, the same form with
    K = diag(mu, 0).  Re(w)/k is the phase speed; -Im(w) the decay rate (1/s).
    """
    k = 2 * math.pi / wavelength
    Mr = np.array([[d['rho'], d['rho_f']], [d['rho_f'], d['m']]])
    Dr = np.array([[0.0, 0.0], [0.0, d['b']]])
    res = {}
    for name, Kp in (('P', np.array([[d['lam_u'] + 2 * d['mu'], d['alpha'] * d['M']],
                                      [d['alpha'] * d['M'], d['M']]])),
                     ('S', np.array([[d['mu'], 0.0], [0.0, 0.0]]))):
        # w^2 Mr X + i w Dr X - k^2 Kp X = 0; with Y = w X: [0 I; Minv k2 K, -i Minv D] [X; Y] = w [X; Y]
        Mi = np.linalg.inv(Mr)
        A = np.block([[np.zeros((2, 2)), np.eye(2)], [Mi @ (k * k * Kp), -1j * Mi @ Dr]])
        ws = np.linalg.eigvals(A)
        ws = sorted([w for w in ws if w.real >= -1e-9 * abs(w) and abs(w) > 0], key=lambda w: -w.real)
        res[name] = [dict(speed=w.real / k, decay=-w.imag) for w in ws]
    return res


def fitted_gamma(z):
    """Eq. (9.4) with the series below z = 1e-2."""
    if z < 1e-2:
        return 0.5 + z / 12.0 - z ** 3 / 720.0 + z ** 5 / 30240.0
    return 1.0 / (-math.expm1(-z)) - 1.0 / z


def gamma_table(zs):
    rows = []
    for z in zs:
        g = fitted_gamma(z)
        R = (1 - (1 - g) * z) / (1 + g * z)                  # eq. (9.3)
        rows.append(dict(z=z, gamma_w=g, R=R, exp_minus_z=math.exp(-z),
                         R_abs_err=abs(R - math.exp(-z))))
    zc = 1e-2
    closed = 1.0 / (-math.expm1(-zc)) - 1.0 / zc
    series = 0.5 + zc / 12.0 - zc ** 3 / 720.0 + zc ** 5 / 30240.0
    return rows, abs(closed - series)


def tpv5_start(d):
    """Section 7.4: sigma_n' at t = 0 and the overstressed patch."""
    m0, w0, x0 = 5e7, 1e-7, np.array([0.0, -1000.0, -7500.0])
    p_hyp = m0 * math.exp(-w0 * 1000.0 ** 2)              # nearest fault point (0, 0, -7.5 km)
    sn, tau0, mus = 120e6, 70e6, 0.677
    p_crit = sn - tau0 / mus
    r3 = math.sqrt(math.log(m0 / p_crit) / w0)
    r_fault = math.sqrt(r3 ** 2 - 1000.0 ** 2)
    # rev 2.2 instant injection: sigma_n' drops by (1 - alpha) p only
    sn_inst = sn - (1 - d['alpha']) * p_hyp
    return dict(p_init_hypocentre=p_hyp, sigma_n_eff_seissol=sn - p_hyp, strength=mus * (sn - p_hyp),
                p_crit=p_crit, overstressed_radius_on_fault=r_fault,
                sigma_n_eff_instant_injection=sn_inst)


def patch_fields(A, B, S0=-1e6, u0=1e-3, V=0.0):
    """PT-1 (V = 0) and PT-2 (uniform V): exact piecewise-linear 1D fields across x = 0.

    u = (u_x(x), 0, 0), w = (w_x(x), 0, 0); material A for x < 0, B for x > 0; drained
    outer faces (p = 0) and exact u on the Dirichlet faces.  Static equilibrium needs
    sigma_xx and p uniform; p = 0 everywhere, so zeta = V on each side.
    """
    out = {}
    for side, d in (('A', A), ('B', B)):
        du = S0 / (d['lam'] + 2 * d['mu'])                 # sigma_xx = (lam + 2 mu) u' when p = 0
        dw = V - d['alpha'] * du                            # zeta = alpha u' + w' = V
        out[side] = dict(du_dx=du, dw_dx=dw,
                         sigma_xx=(d['lam_u'] + 2 * d['mu']) * du + d['alpha'] * d['M'] * (dw - V),
                         p=d['M'] * (V - (d['alpha'] * du + dw)))
    # u(0) = u0 and w continuous at 0 (w(0) = 0); the rev 2.2 penalized jump is [[(alpha u + w).n]]
    out['u_at_interface'] = u0
    out['old_form_jump'] = (A['alpha'] - B['alpha']) * u0
    out['S0'] = S0
    out['V'] = V
    return out


# The material sets of the header (step P1.1) and the scalars written for each: the ten inputs,
# then eq. (4.4), (8.4) and the speeds of (4.6).
HEADER_SETS = ('tpv5_poro', 'tpv205_equivalent', 'vf3_contrast')
HEADER_KEYS = ('mu', 'lam', 'Kf', 'Ks', 'kappa', 'phi', 'rho_f', 'rho_s', 'T', 'eta',
               'Kd', 'alpha', 'M', 'lam_u', 'Ku', 'B', 'rho', 'm', 'b', 'tau', 'f_c',
               'c_fast', 'c_slow', 'c_s_hf', 'c_s_lf', 'c_p_undrained', 'c_p_drained', 'diffusivity')
# The patch-test fields of plan section 6.3 (step P2.2, tests U2-6 and U2-7): per side A (x < 0,
# tpv5-poro) and B (x > 0, TPV205-equivalent), then the shared values.
PATCH_SETS = ('patch_PT1', 'patch_PT2')
PATCH_SIDE_KEYS = ('du_dx', 'dw_dx', 'sigma_xx', 'p')
PATCH_KEYS = ('u_at_interface', 'old_form_jump', 'S0', 'V')


def write_header(refs, path):
    """poro_test_refs.hpp: every value with 17 significant digits (exact round trip)."""
    lines = ['// poro_test_refs.hpp: GENERATED by tests/refs/poro_test_refs.py --header.  Do not edit;',
             '// `make check-poro-refs` regenerates it and fails on any difference.',
             '//',
             '// Reference values for the poroelastic unit tests (plan rev 3, section 6.2), from the',
             '// definitions of poro-derivation-rev1.2, computed with NumPy and SciPy independently of',
             '// the C++ code.  SI units.',
             '',
             '#ifndef MFEM_SEAS_PORO_TEST_REFS_HPP',
             '#define MFEM_SEAS_PORO_TEST_REFS_HPP',
             '',
             'namespace poro_refs',
             '{']
    for name in HEADER_SETS:
        d = refs[name]
        lines.append(f'namespace {name}')
        lines.append('{')
        for k in HEADER_KEYS:
            lines.append(f'constexpr double {k} = {float(d[k]):.16e};')
        lines.append(f'}} // namespace {name}')
    for name in PATCH_SETS:
        d = refs[name]
        lines.append(f'namespace {name}')
        lines.append('{')
        for side in ('A', 'B'):
            for k in PATCH_SIDE_KEYS:
                lines.append(f'constexpr double {side}_{k} = {float(d[side][k]):.16e};')
        for k in PATCH_KEYS:
            lines.append(f'constexpr double {k} = {float(d[k]):.16e};')
        lines.append(f'}} // namespace {name}')
    lines += ['} // namespace poro_refs', '', '#endif // MFEM_SEAS_PORO_TEST_REFS_HPP', '']
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines))


EIGEN_KEYS = ('modes_1km', 'dispersion_tpv5')   # refs.json entries computed by eigensolves


def json_differences(a, b, path=''):
    """(path, committed, generated) for every value of a that differs from b, including a
    difference in keys, list lengths or types."""
    if isinstance(a, dict) and isinstance(b, dict):
        if list(a) != list(b):
            return [(path, list(a), list(b))]
        return [d for k in a for d in json_differences(a[k], b[k], f'{path}/{k}')]
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [(path, len(a), len(b))]
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in json_differences(x, y, f'{path}[{i}]')]
    return [] if a == b and type(a) is type(b) else [(path, a, b)]


def check(refs_dir, generate):
    """Regenerate refs.json and the header into a temporary folder and compare with refs_dir."""
    with tempfile.TemporaryDirectory() as tmp:
        generate(os.path.join(tmp, 'refs.json'), os.path.join(tmp, 'poro_test_refs.hpp'))
        ok = True
        for name, committed in (('refs.json', os.path.join(refs_dir, 'refs', 'refs.json')),
                                ('poro_test_refs.hpp', os.path.join(refs_dir, 'poro_test_refs.hpp'))):
            new = os.path.join(tmp, name)
            if filecmp.cmp(new, committed, shallow=False):
                print(f'  PASSED: check-poro-refs: {name} regenerates byte for byte')
                continue
            if name == 'refs.json':
                diffs = json_differences(json.load(open(committed)), json.load(open(new)))
                top = [p.split('/')[1].split('[')[0] if '/' in p else '' for p, _, _ in diffs]
                eig = [d for d, t in zip(diffs, top) if t in EIGEN_KEYS]
                for p, x, y in diffs:
                    size = (f' (|difference| {abs(x - y):.2e}'
                            + (f', relative {abs(x - y) / abs(x):.2e}' if x != 0 else '') + ')'
                            if all(isinstance(v, float) for v in (x, y)) else '')
                    print(f'  {p}: committed {x!r}, generated {y!r}{size}'
                          + (', eigenvalue result: not checked' if (p, x, y) in eig else ''))
                if diffs and len(eig) == len(diffs):
                    print(f'  PASSED: check-poro-refs: refs.json regenerates with every value equal outside '
                          f'the eigenvalue results; {len(eig)} eigenvalue results differ (listed above)')
                    continue
                ok = False
                print('FAILED: check-poro-refs: refs.json differs from the regenerated file'
                      + (f' in {len(diffs) - len(eig)} values outside the eigenvalue results' if diffs
                         else '; the values agree, so the layout differs'))
                continue
            ok = False
            print(f'FAILED: check-poro-refs: {name} differs from the regenerated file')
            a, b = open(committed).read().splitlines(), open(new).read().splitlines()
            for i, (x, y) in enumerate(zip(a, b)):
                if x != y:
                    print(f'  line {i + 1}: committed: {x}')
                    print(f'  line {i + 1}: generated: {y}')
            if len(a) != len(b):
                print(f'  {len(a)} lines committed, {len(b)} generated')
        return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json')
    ap.add_argument('--header')
    ap.add_argument('--check', metavar='DIR', help='the tests/ folder holding refs/refs.json and the header')
    args = ap.parse_args()
    t5 = derive(6.40e9, 5.28e9, 2.50e9, 2.00e10, 6.00e-13, 0.4, 1040.0, 2080.0, 2.0, 1e-3)
    t205 = tpv205_equivalent()
    vf3 = tpv205_equivalent(phi=0.2, kap=1e-13)
    refs = dict(tpv5_poro=t5, tpv205_equivalent=t205, vf3_contrast=vf3)
    for name, d in refs.items():
        print(f"== {name}")
        for k in ('alpha', 'M', 'lam_u', 'Ku', 'B', 'rho', 'm', 'b', 'tau', 'f_c', 'c_fast', 'c_slow',
                  'c_s_hf', 'c_s_lf', 'c_p_undrained', 'c_p_drained', 'diffusivity'):
            print(f"  {k:14s} {d[k]:.12g}")
        print(f"  tau vs eig rel {abs(d['tau'] - d['tau_from_eig']) / d['tau']:.2e}; "
              f"M forms rel {d['M_textbook_rel_diff']:.2e}")
        for h, p in ((1000, 1), (200, 1), (200, 2), (200, 3), (100, 2), (100, 3)):
            dt = 0.5 * h / (p * p * d['c_fast'])
            print(f"  h={h:5d} p={p}: dt(cfl 0.5) = {dt:.3e} s, z = {dt / d['tau']:.3g}, "
                  f"gamma_w = {fitted_gamma(dt / d['tau']):.6f}")
    print(f"  TPV205 set: lam_u rel err vs TPV205 lambda {t205['lam_u_rel_err']:.2e}")
    print(f"  V-F3 set:   lam_u rel err vs TPV205 lambda {vf3['lam_u_rel_err']:.2e}")
    disp = dispersion(t5, [1.0, 10.0, t5['f_c'], 100 * t5['f_c']])
    print("== tpv5-poro dispersion (speed m/s, Q)")
    for r in disp:
        print(f"  f={r['f']:.4g} Hz: fastP {r['fastP_speed']:.2f} (Q {r['fastP_Q']:.3g}), "
              f"slowP {r['slowP_speed']:.4g} (Q {r['slowP_Q']:.3g}), S {r['S_speed']:.2f} (Q {r['S_Q']:.3g})")
    rows, junction = gamma_table([1e-8, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 2.6e3, 2.1e5, 1e8])
    print("== fitted gamma_w")
    for r in rows:
        print(f"  z={r['z']:.3g}: gamma_w={r['gamma_w']:.15f}  |R - e^-z| = {r['R_abs_err']:.2e}")
    print(f"  series vs closed form at z = 1e-2: {junction:.2e}")
    st = tpv5_start(t5)
    print("== tpv5-poro start")
    for k, v in st.items():
        print(f"  {k:32s} {v:.6g}")
    pt1 = patch_fields(t5, t205)
    pt2 = patch_fields(t5, t205, V=1e-4)
    print("== patch tests (A = tpv5-poro for x < 0, B = TPV205-equivalent for x > 0)")
    for name, pt in (('PT-1', pt1), ('PT-2', pt2)):
        print(f"  {name}: V = {pt['V']:g}; A: u' {pt['A']['du_dx']:.6e}, w' {pt['A']['dw_dx']:.6e}, "
              f"sxx {pt['A']['sigma_xx']:.6e}, p {pt['A']['p']:.2e}; B: u' {pt['B']['du_dx']:.6e}, "
              f"w' {pt['B']['dw_dx']:.6e}, sxx {pt['B']['sigma_xx']:.6e}, p {pt['B']['p']:.2e}; "
              f"rev2.2 jump [[(alpha u + w).n]] = {pt['old_form_jump']:.4e} m")
    print("== periodic-box modes, wavelength 1 km (speed m/s, decay 1/s)")
    modes = {}
    for name, d in (('tpv5-poro', t5), ('tpv5-poro kappa=1e-7', derive(6.40e9, 5.28e9, 2.50e9, 2.00e10,
                    1e-7, 0.4, 1040.0, 2080.0, 2.0, 1e-3))):
        md = modes_real_k(d, 1000.0)
        modes[name] = md
        for kind in ('P', 'S'):
            print(f"  {name:22s} {kind}: " + ', '.join(f"{r['speed']:.2f}/{r['decay']:.4g}" for r in md[kind]))
    refs.update(modes_1km=modes, dispersion_tpv5=disp, gamma=rows, gamma_junction=junction, tpv5_start=st,
                patch_PT1=pt1, patch_PT2=pt2)

    def generate(json_path, header_path):
        if json_path:
            with open(json_path, 'w') as fh:
                json.dump(refs, fh, indent=1, default=float)
        if header_path:
            write_header(refs, header_path)

    if args.check:
        sys.exit(0 if check(args.check, generate) else 1)
    generate(args.json, args.header)


if __name__ == '__main__':
    main()
