#!/bin/bash
# Mirror-symmetric TPV26/TPV27 meshes (rheology plan rev 1, phase R1; Chunhui Zhao,
# 2026-10-20): gmsh meshes the y >= 0 half (tpv26/mesh/tpv2627_half.geo) and
# scripts/mirror_mesh.py reflects it across Y = 0, as scripts/make_mirror_meshes.sh does for
# TPV104/TPV205.  Run from miniapps/seas with gmsh and python3 on PATH (fault element sizes
# in m; extra gmsh -setnumber options, e.g. box_halfwidth, go in $GMSH_OPTS):
#   bash tpv26/mesh/make_tpv2627_meshes.sh 1000 200
# Writes tpv26/mesh/tpv2627_mirror_<res>m.msh (gmsh and mirror logs next to it) and prints
# "sha256  path" for each mesh, also to $SHA_FILE when it is set.  gmsh output depends on the
# platform, so the sha256 identifies the mesh used.
set -eo pipefail
cd "$(dirname "$0")/../.."
command -v gmsh > /dev/null || { echo "ERROR: gmsh not on PATH"; exit 2; }
[ $# -gt 0 ] || { echo "usage: bash tpv26/mesh/make_tpv2627_meshes.sh <res_m> [<res_m> ...]"; exit 2; }
echo "gmsh $(gmsh --version 2>&1), $(python3 --version 2>&1)"
sha() { if command -v sha256sum > /dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
meshes=()
for res in "$@"; do
   base=tpv26/mesh/tpv2627_mirror_${res}m
   t=$(date +%s)
   # shellcheck disable=SC2086
   gmsh -3 tpv26/mesh/tpv2627_half.geo -setnumber lc_fault "$res" $GMSH_OPTS -format msh22 \
        -o "$base.half.msh" > "$base.gmsh.log" 2>&1 \
      || { echo "ERROR: gmsh failed at $res m (log $base.gmsh.log)"; tail -5 "$base.gmsh.log"; exit 1; }
   # gmsh exits 0 when a .geo command fails (a station point it cannot embed), so read its log too.
   if grep -q "^Error" "$base.gmsh.log"; then
      echo "ERROR: gmsh reported errors at $res m (log $base.gmsh.log)"; grep "^Error" "$base.gmsh.log"; exit 1
   fi
   python3 scripts/mirror_mesh.py "$base.half.msh" "$base.msh" --fault-attr 103 | tee "$base.mirror.log"
   rm -f "$base.half.msh"
   echo "$base.msh: $(( $(date +%s) - t )) s"
   meshes+=("$base.msh")
done
if [ -n "$SHA_FILE" ]; then sha "${meshes[@]}" | tee "$SHA_FILE"; else sha "${meshes[@]}"; fi
