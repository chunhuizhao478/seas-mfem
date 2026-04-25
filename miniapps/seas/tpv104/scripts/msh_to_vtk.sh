#!/usr/bin/env bash
# Convert one or more Gmsh .msh files to legacy VTK for ParaView viewing.
#
# Usage:
#   ./msh_to_vtk.sh file1.msh [file2.msh ...]
#   ./msh_to_vtk.sh ../mesh/*.msh
#
# Output: file.vtk next to each file.msh (binary VTK legacy format).
#
# In ParaView:
#   1. File -> Open <file>.vtk
#   2. The cell data array "CellEntityIds" carries Gmsh physical group
#      IDs.  Use the Threshold filter:
#        CellEntityIds == 1 -> free surface (Z=0)
#        CellEntityIds == 3 -> fault face (Y=0 interior)
#        CellEntityIds == 5 -> absorbing outer boundaries
#   3. For a symmetry-check viewport, slice along Y=0 and color by
#      mesh quality / cell volume to spot non-mirror tet pairs.
#
# Requires: gmsh (in conda env `pythonenv` per CLAUDE.md).

set -euo pipefail

if [ "$#" -lt 1 ]; then
   echo "Usage: $0 <file1.msh> [file2.msh ...]" >&2
   exit 1
fi

if ! command -v gmsh >/dev/null 2>&1; then
   echo "ERROR: gmsh not found.  Activate the pythonenv conda env first:" >&2
   echo "         conda activate pythonenv" >&2
   exit 1
fi

for msh in "$@"; do
   if [ ! -s "$msh" ]; then
      echo "  skip (empty/missing): $msh" >&2
      continue
   fi
   vtk="${msh%.msh}.vtk"
   echo "Converting $msh -> $vtk ..."
   gmsh "$msh" -o "$vtk" -3 -bin 2>&1 | grep -E "Writing|Done|Error" || true
done
echo "Done.  Open .vtk files in ParaView."
