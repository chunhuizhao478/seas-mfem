#!/bin/zsh
# ParaView views of the two gate-closed products.
#
# Three per mesh, cheapest first (mesh_xdmf.py):
#   _fault    dynamic-rupture triangles only, de-duplicated (a PUML stores each
#             fault triangle twice, once per side). ~2 % of the cells and the
#             thing you actually want to look at.
#   _surface  fault + free surface + absorbing hull, carrying a `bc` scalar, so
#             the domain box and the lid show up with the fault.
#   _full     a ~2 KB XDMF that POINTS AT the .puml.h5 -- zero bytes copied, but
#             ParaView then loads the whole mesh. Fine for the 38.8M
#             intermediate, heavy going for the 157.8M heavy: open _surface or
#             _fault there unless you have the RAM.
#
# The XDMF resolves <DataItem> paths RELATIVE TO ITSELF, so each .xdmf must stay
# next to the .h5 it names; _full's relative path back to ../results is baked in
# by mesh_xdmf.py, so do not move these files independently.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
X=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/.claude/worktrees/shakeoutbox-gate-close/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_small/code/mesh_xdmf.py
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
R=$A/meshing_shakeoutbox_gate/results
V=$A/meshing_shakeoutbox_gate/view
mkdir -p $V
cd $V

INT=$R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5
HVY=$R/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5

for m in fault surface; do
  echo "########## intermediate $m"
  $PY $X --mesh $INT --mode $m --out $V/intermediate_$m.xdmf
done
echo "########## intermediate full"
$PY $X --mesh $INT --mode wrap --out $V/intermediate_full.xdmf

for m in fault surface; do
  echo "########## heavy $m"
  $PY $X --mesh $HVY --mode $m --out $V/heavy_$m.xdmf
done
echo "########## heavy full"
$PY $X --mesh $HVY --mode wrap --out $V/heavy_full.xdmf
echo ALL DONE
