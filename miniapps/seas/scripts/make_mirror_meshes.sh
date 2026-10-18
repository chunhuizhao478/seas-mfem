#!/bin/bash
# Mirror-symmetric meshes for the Nitsche TPV104 / TPV205 runs (Chunhui Zhao, 2026-10-17):
# gmsh meshes the y >= 0 half (tpv104/mesh/tpv104_half.geo with the 22 km fault and its
# velocity-strengthening border, tpv205/mesh/tpv2053d_half.geo) and scripts/mirror_mesh.py
# reflects it across Y = 0, so the tets on the two sides of every fault face are mirror images.
# Run from miniapps/seas with gmsh and python3 on PATH (resolutions in m):
#   bash scripts/make_mirror_meshes.sh 200 100
# Writes tpv104/mesh/tpv104_mirror_ext22_<res>m.msh and tpv205/mesh/tpv2053d_mirror_<res>m.msh
# (gmsh and mirror logs next to them) and prints "sha256  path" for each mesh, also to $SHA_FILE
# when it is set.  gmsh output depends on the platform, so the sha256 identifies the mesh used.
set -eo pipefail
cd "$(dirname "$0")/.."
command -v gmsh > /dev/null || { echo "ERROR: gmsh not on PATH"; exit 2; }
echo "gmsh $(gmsh --version 2>&1), $(python3 --version 2>&1)"
sha() { if command -v sha256sum > /dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
meshes=()
for res in "$@"; do
   for spec in "tpv104/mesh/tpv104_half.geo tpv104/mesh/tpv104_mirror_ext22_${res}m 3" \
               "tpv205/mesh/tpv2053d_half.geo tpv205/mesh/tpv2053d_mirror_${res}m 103"; do
      read -r geo base fault_attr <<< "$spec"
      t=$(date +%s)
      gmsh -3 "$geo" -setnumber lc_fault "$res" -format msh22 -o "$base.half.msh" > "$base.gmsh.log" 2>&1 \
         || { echo "ERROR: gmsh failed on $geo at $res m (log $base.gmsh.log)"; tail -5 "$base.gmsh.log"; exit 1; }
      python3 scripts/mirror_mesh.py "$base.half.msh" "$base.msh" --fault-attr "$fault_attr" \
         | tee "$base.mirror.log"
      rm -f "$base.half.msh"
      echo "$base.msh: $(( $(date +%s) - t )) s"
      meshes+=("$base.msh")
   done
done
if [ -n "$SHA_FILE" ]; then sha "${meshes[@]}" | tee "$SHA_FILE"; else sha "${meshes[@]}"; fi
