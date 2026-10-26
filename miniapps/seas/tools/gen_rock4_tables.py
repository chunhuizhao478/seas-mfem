#!/usr/bin/env python3
# Usage: python3 -I gen_rock4_tables.py rkc_tableaus_rock4.jl solver/rock4_tables.hpp LICENSE.md
# Writes solver/rock4_tables.hpp from OrdinaryDiffEq.jl's rkc_tableaus_rock4.jl
# (numbers copied as text, so no value changes in the conversion).
import re, sys, hashlib
path, out = sys.argv[1], sys.argv[2]
raw = open(path, "rb").read()
src = raw.decode()
sha = hashlib.sha256(raw).hexdigest()
num = re.compile(r"[-+]?\d+\.\d*(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?")
def section(name, nxt):
    i = src.index("    %s = " % name)
    j = src.index("    %s = " % nxt, i) if nxt else src.index("    _fpa = ", i)
    return num.findall(src[i:j].split("=", 1)[1])
ms = section("ms", "fpa"); fpa = section("fpa", "fpb"); fpb = section("fpb", "fpbe")
fpbe = section("fpbe", "recf"); recf = section("recf", None)
assert len(ms) == 50 and len(fpa) == 300 and len(fpb) == 200 and len(fpbe) == 250
assert len(recf) == sum(2 * int(m) - 1 for m in ms) == 4382
def fmt(x):
    return x if ("." in x or "e" in x or "E" in x) else x + ".0"
def rows(vals, per, indent="   "):
    lines = []
    for k in range(0, len(vals), per):
        lines.append(indent + ", ".join(fmt(v) for v in vals[k:k + per]) + ",")
    return "\n".join(lines)
mit = open(sys.argv[3]).read().strip().splitlines()
with open(out, "w") as f:
    f.write("""// rock4_tables.hpp — coefficients of ROCK4 (A. Abdulle, "Fourth order
// Chebyshev methods with recurrence relation", SIAM J. Sci. Comput. 23(6),
// 2041-2054, 2002), for solver/ts_rock4.hpp.  Copied as text by
// tools/gen_rock4_tables.py, unchanged, from OrdinaryDiffEq.jl:
//   https://raw.githubusercontent.com/SciML/OrdinaryDiffEq.jl/master/
//   lib/OrdinaryDiffEqStabilizedRK/src/rkc_tableaus_rock4.jl
//   fetched 2026-10-25, sha256 %s
// whose tables are those of Abdulle's rock4.f.  That file's license:
//
""" % sha)
    for line in mit:
        f.write(("// " + line).rstrip() + "\n")
    f.write("""//
// Layout (as rock4.f and the Julia file):
//   kRock4Ms[i]       degree s_i of the Chebyshev part of method i (50
//                     methods, s_i = 1..148); the method has s_i + 4 stages.
//   kRock4Fpa[i][6]   finishing procedure a21, a31, a32, a41, a42, a43
//   kRock4Fpb[i][4]   its weights b1..b4 (order 4)
//   kRock4Fpbe[i][5]  the embedded weights (order 3), the fifth on f(y_n+1)
//   kRock4Recf[]      for each method in turn, 2 s_i - 1 values: mu_1, then
//                     (mu_j, kappa_j) for j = 2..s_i

#ifndef MFEM_SEAS_ROCK4_TABLES_HPP
#define MFEM_SEAS_ROCK4_TABLES_HPP

namespace mfem
{
namespace seas
{

constexpr int kRock4Methods = 50;
constexpr int kRock4RecfSize = 4382;

constexpr int kRock4Ms[kRock4Methods] =
{
%s
};

constexpr double kRock4Fpa[kRock4Methods][6] =
{
%s
};

constexpr double kRock4Fpb[kRock4Methods][4] =
{
%s
};

constexpr double kRock4Fpbe[kRock4Methods][5] =
{
%s
};

constexpr double kRock4Recf[kRock4RecfSize] =
{
%s
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_ROCK4_TABLES_HPP
""" % ("\n".join("   " + ", ".join(ms[k:k + 10]) + "," for k in range(0, 50, 10)),
       "\n".join("   {" + ", ".join(fmt(v) for v in fpa[6*k:6*k+6]) + "}," for k in range(50)),
       "\n".join("   {" + ", ".join(fmt(v) for v in fpb[4*k:4*k+4]) + "}," for k in range(50)),
       "\n".join("   {" + ", ".join(fmt(v) for v in fpbe[5*k:5*k+5]) + "}," for k in range(50)),
       rows(recf, 3)))
print("wrote", out, sha)
