#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/nitsche_finite_width_decks.sh — test-2d5-dyn-driver: the
# [finite_width] table of seas_nitsche_dyn_driver (the 2.5D note's Level 2
# test, elongated2p5d/; Chunhui Zhao, 2026-10-30).  Copies of the TPV205-2D
# deck (tpv205/configs/tpv205_2d_nitsche.toml, its committed 100 m mesh)
# with a [finite_width] table appended:
#   - each refusal must exit non-zero with its message (zero W, calibration
#     and c together, an unknown calibration, an unknown key, LTS, the
#     energy balance; a 3-D mesh when MESHES holds tpv2053d_1000m.msh, as
#     tests/2d_inplane/parity_3d.sh);
#   - the three accepted tables (wa2020, lehner, c as a number) must run
#     0.02 s and report W, c, beta and the step's cut, with
#     beta = 4 c^2 mu / W^2 to 6 digits;
#   - the run with the term must change the station output against the run
#     without it (the term reaches the update).
# The deck without the table is checked bit for bit against the binary
# before the option by tests/2d_inplane/parity_3d.sh (3-D) and by the
# review's 2-D run, not here.
#
# Usage (from miniapps/seas, after building seas_nitsche_dyn_driver):
#   MPIEXEC=mpirun MPIEXEC_NP=-np tests/2d_inplane/nitsche_finite_width_decks.sh
# Exit code 0 = PASSED.

set -u
export LC_ALL=C
MPI=${MPIEXEC:-mpirun}
NP=${MPIEXEC_NP:--np}
BASE=tpv205/configs/tpv205_2d_nitsche.toml
BASE3=tpv205/configs/tpv205_nitsche.toml
T=${TMPDIR:-/tmp}/seas_test_2d5_dyn
rm -rf "$T"
mkdir -p "$T"
fail=0

# deck NAME BASE TEXT: BASE with TEXT (printf %b) appended
deck() { cp "$2" "$T/$1.toml" && printf '%b' "$3" >> "$T/$1.toml"; }
run() { $MPI $NP 1 ./seas_nitsche_dyn_driver --config "$T/$1.toml" --output-dir "$T/$1_out" \
            --tfinal 0.02 --log-every 1000 "${@:2}" 2>&1; }

# refuse NAME MESSAGE [driver options]: the deck $T/NAME.toml must be refused
refuse()
{
   out=$(run "$1" "${@:3}")
   rc=$?
   if [ $rc -ne 0 ] && printf '%s' "$out" | grep -qF -- "$2"; then
      echo "  PASSED: refusal $1: exit $rc, \"$2\""
   else
      printf '%s\n' "$out" | tail -15
      echo "  FAILED: refusal $1: exit $rc, \"$2\" not found"
      fail=1
   fi
}

H='\n[finite_width]\n'
deck zero "$BASE" "${H}W_m = 0.0\n"
refuse zero "[finite_width].W_m = 0.0000000000000000e+00 must be > 0 m"
deck both "$BASE" "${H}W_m = 15.0e3\ncalibration = \"wa2020\"\nc = 1.2\n"
refuse both "[finite_width] takes calibration (a name) or c (a number), not both"
deck name "$BASE" "${H}W_m = 15.0e3\ncalibration = \"surface\"\n"
refuse name "[finite_width].calibration = \"surface\"; use \"wa2020\""
deck key "$BASE" "${H}W_m = 15.0e3\nwidth = 1.0\n"
refuse key "[finite_width]: unknown key 'width' (W_m, calibration, c)"
deck lts "$BASE" "${H}W_m = 15.0e3\n"
refuse lts "[finite_width] with LTS: the LTS run integrates its own element matrices" --lts rate2
awk '{ sub(/^energy_balance = false/, "energy_balance = true") } 1' "$BASE" > "$T/energy.toml"
grep -q '^energy_balance = true' "$T/energy.toml" || { echo "  FAILED: energy deck not made"; fail=1; }
printf '%b' "${H}W_m = 15.0e3\n" >> "$T/energy.toml"
refuse energy "[finite_width] with [diagnostics] energy_balance: the energy terms carry no beta"
if [ -n "${MESHES:-}" ] && [ -f "$MESHES/tpv2053d_1000m.msh" ]; then
   awk '{ sub(/^paraview_fault *= *"hdf5"/, "paraview_fault = \"off\"") } 1' "$BASE3" > "$T/three.toml"
   printf '%b' "${H}W_m = 15.0e3\n" >> "$T/three.toml"
   refuse three "[finite_width]: the 2.5D term is for 2-D map-view meshes; this mesh is 3-D" \
          --mesh "$MESHES/tpv2053d_1000m.msh"
else
   echo "  SKIPPED: refusal three (3-D mesh): set MESHES to the folder with tpv2053d_1000m.msh"
fi

# accept NAME TEXT C NAME_IN_LOG: runs 0.02 s; the log's c and beta
# (TPV205-2D's [material_constant_fallback] mu)
mu=32038120320.0
accept()
{
   deck "$1" "$BASE" "$2"
   out=$(run "$1")
   rc=$?
   line=$(printf '%s' "$out" | grep -F '2.5D finite width ([finite_width]): W = 15000 m' | head -1)
   want_c=$(awk -v c="$3" 'BEGIN { printf "%.6g", c }')
   want_b=$(awk -v c="$3" -v mu="$mu" 'BEGIN { printf "%.6g", 4 * c * c * mu / 15000 / 15000 }')
   if [ $rc -eq 0 ] && printf '%s' "$line" | grep -qF "c = $want_c ($4), beta = $want_b Pa/m^2" \
      && printf '%s' "$line" | grep -qF "the CFL step times (1 - "; then
      echo "  PASSED: accepted $1 (0.02 s): $(printf '%s' "$line" | sed 's/^ *//')"
   else
      printf '%s\n' "$out" | tail -20
      echo "  FAILED: accepted $1: exit $rc; want c = $want_c ($4), beta = $want_b"
      fail=1
   fi
}
# wa2020: c = (1 - nu) pi / 2 with nu from the deck's lambda and mu (TPV205:
# nu = 0.250022, c = 1.17806)
lam=32043759360.0
c_wa=$(awk -v l="$lam" -v m="$mu" 'BEGIN { printf "%.17g", (1 - l / (2 * (l + m))) * atan2(0, -1) / 2 }')
accept wa2020 "${H}W_m = 15.0e3\n" "$c_wa" wa2020
accept lehner "${H}W_m = 15.0e3\ncalibration = \"lehner\"\n" 1.2732395447351628 lehner
accept number "${H}W_m = 15.0e3\nc = 1.5707963267948966\n" 1.5707963267948966 deck

# the term reaches the update: the stations differ from the run without it
cp "$BASE" "$T/plain.toml"
out=$(run plain)
rc=$?
if [ $rc -ne 0 ]; then
   printf '%s\n' "$out" | tail -15
   echo "  FAILED: the run without [finite_width]: exit $rc"
   fail=1
else
   n=0; same=0
   for f in "$T"/plain_out/*station*; do
      g="$T/wa2020_out/$(basename "$f")"
      [ -f "$g" ] || continue
      n=$((n + 1))
      cmp -s <(grep -v '^#' "$f") <(grep -v '^#' "$g") && same=$((same + 1))
   done
   if [ $n -gt 0 ] && [ $same -lt $n ]; then
      echo "  PASSED: the term changes the stations ($((n - same)) of $n files differ from the run without it)"
   else
      echo "  FAILED: the term changes the stations: $n files compared, $same identical"
      fail=1
   fi
fi

rm -rf "$T"
exit $fail
