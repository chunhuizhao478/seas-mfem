#!/usr/bin/env python3
"""Compare the poro verification runs with their closed forms.

  python3 check_poro_verification.py undrained <undrained_column_c1.csv>
  python3 check_poro_verification.py terzaghi  <terzaghi_column_c1.csv>
  python3 check_poro_verification.py deleeuw   <deleeuw_cylinder_c1.csv>

Parameters are those of undrained_column.toml / terzaghi_column.toml.
"""
import csv, math, sys

E, nu, rho = 10e9, 0.3, 2600.0
K_s, K_f, phi0 = 50e9, 2.24e9, 0.2
L, W, P = 0.01, 0.0005, 1e6
K = E / 3 / (1 - 2 * nu); G = E / 2 / (1 + nu)
alpha = 1 - K / K_s
M = 1 / (phi0 / K_f + (alpha - phi0) / K_s)
H = K + 4 * G / 3; Hu = H + alpha**2 * M
cu = math.sqrt(Hu / rho)

kind, path = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(path)))
t = [float(r["t"]) for r in rows]
col = lambda k: [float(r[k]) for r in rows]
fail = 0
def check(name, got, exp, rtol):
    global fail
    ok = abs(got - exp) <= rtol * abs(exp)
    fail += not ok
    print("  %s %s: %.6g, expected %.6g (%.2f %%, tolerance %.1f %%)" %
          ("PASS" if ok else "FAIL", name, got, exp, 100 * (got / exp - 1), 100 * rtol))

if kind == "undrained":
    print("undrained column: H_u = %.4g Pa (drained %.4g), c_u = %.1f m/s, 2L/c_u = %.4g s" % (Hu, H, cu, 2 * L / cu))
    wc = col("confinement_work")
    i = max(range(len(t)), key=lambda k: wc[k] if t[k] < 3 * L / cu else -1)
    check("max confinement work 2 P^2 L W / H_u", wc[i], 2 * P * P * L * W / Hu, 0.01)
    check("time of the maximum 2L/c_u", t[i], 2 * L / cu, 0.02)
    ke = col("kinetic")
    j = min((k for k in range(len(t)) if 1.5 * L / cu < t[k] < 2.5 * L / cu), key=lambda k: ke[k])
    check("time of the kinetic energy minimum 2L/c_u", t[j], 2 * L / cu, 0.02)
    # at t = L/c_u the whole column moves at v = P/(rho c_u): KE = 1/2 rho v^2 L W
    k1 = min(range(len(t)), key=lambda k: abs(t[k] - L / cu))
    check("kinetic energy at L/c_u, P^2 L W / (2 H_u)", ke[k1], P * P * L * W / (2 * Hu), 0.03)
    err = max(abs(e) for e in col("energy_error"))
    print("  max |energy_error| / max work = %.2e" % (err / max(wc)))
    fail += err / max(wc) > 1e-3
elif kind == "deleeuw":
    # p(0, s) = alpha P (1 - 1/I0(x)) / (H S s D(x)), D = 1 - m + m e - 2 m e I1(x)/(x I0(x)),
    # x = a sqrt(s/c). Poles: s = 0 (zero residue) and s_n = -c y_n^2 / a^2 with
    # f(y) = (1 - m + m e) y J0(y) - 2 m e J1(y) = 0; residue series below.
    a = 0.01
    k_mu = 3.84e-15 / 1e-3
    S = 1 / M + alpha**2 / H
    c = k_mu / S
    m = G / H
    e = alpha**2 / (H * S)
    def J(n, y, N=2000):   # Bessel's integral, trapezoid (periodic integrand)
        return sum(math.cos(n * math.pi * (k + 0.5) / N - y * math.sin(math.pi * (k + 0.5) / N)) for k in range(N)) / N
    f = lambda y: (1 - m + m * e) * y * J(0, y) - 2 * m * e * J(1, y)
    roots, y0, dy = [], 1e-3, 0.02
    while len(roots) < 60:
        y1 = y0 + dy
        if f(y0) * f(y1) < 0:
            lo, hi = y0, y1
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
            roots.append(0.5 * (lo + hi))
        y0 = y1
    D = lambda y: 1 - m + m * e - 2 * m * e * J(1, y) / (y * J(0, y))
    terms = []
    for y in roots:
        h = 1e-6 * y
        dDdy = (D(y + h) - D(y - h)) / (2 * h)
        dDds = dDdy / (-2 * c * y / a**2)
        sn = -c * y * y / a**2
        terms.append((sn, alpha * P / (H * S) * (1 - 1 / J(0, y)) / (sn * dDds)))
    pc = lambda tt: sum(A * math.exp(sn * tt) for sn, A in terms)
    p_und = alpha * M * P / (Hu - G)
    print("de Leeuw cylinder: c = %.5g m^2/s, undrained centre pressure alpha M P/(H_u - G) = %.6g Pa" % (c, p_und))
    print("  series at T = 1e-4: %.6g Pa (60 roots; first roots %s)" % (pc(1e-4 * a * a / c), ", ".join("%.4f" % y for y in roots[:3])))
    pm = col("p_max")
    tpk = max((0.002 * i * a * a / c for i in range(1, 200)), key=pc)
    print("  Mandel-Cryer peak of the series: %.6g Pa (%.2f %% above undrained) at T = %.3f" % (pc(tpk), 100 * (pc(tpk) / p_und - 1), c * tpk / a**2))
    kpk = max(range(len(t)), key=lambda i: pm[i] if t[i] > 0.01 * a * a / c else 0)
    check("peak centre pressure", pm[kpk], pc(tpk), 0.01)
    check("overshoot above undrained (peak - p_und)", pm[kpk] - p_und, pc(tpk) - p_und, 0.10)
    for Tv in (0.02, 0.05, 0.1, 0.15, 0.2, 0.3):
        tt = Tv * a * a / c
        k = min(range(len(t)), key=lambda i: abs(t[i] - tt))
        check("centre pressure, T = c t / a^2 = %.2f" % Tv, pm[k], pc(t[k]), 0.01)
    wc, err = col("confinement_work"), col("energy_error")
    print("  energy_error at the end / confinement work = %.2e" % (err[-1] / wc[-1]))
else:
    k_mu = 3.84e-15 / 1e-3
    cv = k_mu * M * H / Hu
    p0 = alpha * M * P / Hu
    print("Terzaghi: p0 = %.5g Pa, c_v = %.5g m^2/s, T_v(end) = %.3f" % (p0, cv, cv * t[-1] / L**2))
    pm = col("p_max")
    def series(tt):
        s = 0.0
        for m in range(400):
            a = (2 * m + 1) * math.pi / 2
            s += 2 / a * math.sin(a) * math.exp(-a * a * cv * tt / L**2)
        return p0 * s
    for Tv in (0.05, 0.1, 0.2, 0.3, 0.5):
        tt = Tv * L**2 / cv
        k = min(range(len(t)), key=lambda i: abs(t[i] - tt))
        check("p at the sealed bottom, T_v = %.2f" % Tv, pm[k], series(t[k]), 0.01)
    # settlement: confinement work = P W u_top; u_top from P L/H_u (undrained) to P L/H (drained)
    wc = col("confinement_work")
    def settle(tt):
        s = 0.0
        for m in range(400):
            a = (2 * m + 1) * math.pi / 2
            s += 2 / (a * a) * math.exp(-a * a * cv * tt / L**2)
        return P * L / Hu + (P * L / H - P * L / Hu) * (1 - s)
    for Tv in (0.1, 0.3, 0.5):
        tt = Tv * L**2 / cv
        k = min(range(len(t)), key=lambda i: abs(t[i] - tt))
        check("top settlement (confinement work / P W), T_v = %.2f" % Tv, wc[k] / (P * W), settle(t[k]), 0.01)
    err = col("energy_error")
    print("  energy_error at the end / confinement work = %.2e" % (err[-1] / wc[-1]))
print("FAILED" if fail else "all checks passed")
sys.exit(1 if fail else 0)
