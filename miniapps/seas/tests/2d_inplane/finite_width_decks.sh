#!/bin/bash
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/finite_width_decks.sh — the deck side of test-2d5-width
# ([cycle.finite_width], cycle/finite_width.hpp; Chunhui Zhao, 2026-10-30):
# copies of the map-view deck (config/fd_cycle/mapview_rs_tensor_fd.toml)
# with a [cycle.finite_width] table appended, run through
# seas_fd_cycle_driver.  Each refusal must exit non-zero with its message;
# the two accepted tables (wa2020 + loading_plane, lehner + locked_reference)
# must run 3 QD steps and report W, c, the decay length and the substrate.
#
# Usage (from miniapps/seas, after building seas_fd_cycle_driver):
#   MPIEXEC=mpirun MPIEXEC_NP=-np tests/2d_inplane/finite_width_decks.sh
# Exit code 0 = PASSED.

set -u
MPI=${MPIEXEC:-mpirun}
NP=${MPIEXEC_NP:--np}
MAP=config/fd_cycle/mapview_rs_tensor_fd.toml
BP3=config/fd_cycle/bp3_fd_d60_thrust.toml
T=${TMPDIR:-/tmp}/seas_test_2d5_width
rm -rf "$T"
mkdir -p "$T"
fail=0

# deck NAME BASE TEXT: BASE with TEXT (printf %b) appended
deck() { cp "$2" "$T/$1.toml" && printf '%b' "$3" >> "$T/$1.toml"; }
run() { $MPI $NP 1 ./seas_fd_cycle_driver -d "$T/$1.toml" -o "$T/$1_out" "${@:2}" 2>&1; }

# refuse NAME BASE TEXT MESSAGE
refuse()
{
   deck "$1" "$2" "$3"
   out=$(run "$1")
   rc=$?
   if [ $rc -ne 0 ] && printf '%s' "$out" | grep -qF -- "$4"; then
      echo "  PASSED: deck refusal $1: exit $rc, \"$4\""
   else
      printf '%s\n' "$out" | tail -15
      echo "  FAILED: deck refusal $1: exit $rc, \"$4\" not found"
      fail=1
   fi
}

H='\n[cycle.finite_width]\n'
refuse zero   "$MAP" "${H}W_m = 0.0\nsubstrate = \"loading_plane\"\n" \
       "[cycle.finite_width].W_m = 0.0000000000000000e+00 must be > 0 m"
refuse both   "$MAP" "${H}W_m = 15.0e3\ncalibration = \"wa2020\"\nc = 1.2\nsubstrate = \"loading_plane\"\n" \
       "[cycle.finite_width] takes calibration (a name) or c (a number), not both"
refuse name   "$MAP" "${H}W_m = 15.0e3\ncalibration = \"surface\"\nsubstrate = \"loading_plane\"\n" \
       "[cycle.finite_width].calibration = \"surface\"; use \"wa2020\""
refuse nosub  "$MAP" "${H}W_m = 15.0e3\n" \
       "[cycle.finite_width].substrate is required"
refuse badsub "$MAP" "${H}W_m = 15.0e3\nsubstrate = \"plate\"\n" \
       "use \"loading_plane\" or \"locked_reference\""
refuse key    "$MAP" "${H}W_m = 15.0e3\nwidth = 1.0\nsubstrate = \"loading_plane\"\n" \
       "[cycle.finite_width] has an unknown key 'width'"
refuse lts    "$MAP" "${H}W_m = 15.0e3\nsubstrate = \"loading_plane\"\n\n[cycle.lts]\nmode = \"rate2\"\n" \
       "[cycle.finite_width] with [cycle.lts] mode = \"rate2\""
refuse bp3    "$BP3" "${H}W_m = 15.0e3\nsubstrate = \"locked_reference\"\n" \
       "the 2.5D term is for map-view decks"

# accept NAME TEXT EXPECTED...
accept()
{
   deck "$1" "$MAP" "$2"
   out=$(run "$1" --max-steps 3)
   rc=$?
   ok=1
   [ $rc -eq 0 ] || ok=0
   for want in "${@:3}"; do printf '%s' "$out" | grep -qF -- "$want" || { ok=0; echo "    missing: $want"; }; done
   if [ $ok -eq 1 ]; then
      echo "  PASSED: deck $1 accepted (3 QD steps): $(printf '%s' "$out" | grep -F '2.5D finite width (' | sed 's/^ *//')"
   else
      printf '%s\n' "$out" | tail -20
      echo "  FAILED: deck $1: exit $rc"
      fail=1
   fi
}
accept plane "${H}W_m = 15.0e3\nsubstrate = \"loading_plane\"\n" \
       "2.5D finite width ([cycle.finite_width]): W = 15000 m, c = 1.1781 (wa2020)" \
       "substrate \"loading_plane\" (a locked straight fault on the plane is loaded at"
accept locked "${H}W_m = 15.0e3\ncalibration = \"lehner\"\nsubstrate = \"locked_reference\"\n" \
       "2.5D finite width ([cycle.finite_width]): W = 15000 m, c = 1.27324 (lehner)" \
       "substrate \"locked_reference\" (the loading of the run without the term"

rm -rf "$T"
exit $fail
