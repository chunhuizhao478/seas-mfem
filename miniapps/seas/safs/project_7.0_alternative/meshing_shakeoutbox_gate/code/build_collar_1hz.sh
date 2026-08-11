#!/bin/zsh
# Base collar for the 1 Hz @ p5 heavy parent.
#
# The heavy mesh currently welds its collar onto `safalt_fb200_deep40km_refine2`
# (122,162,105 tets), which fails the 1 Hz gate on 274,299 of its own cells --
# 75 % of that product's total failures.  `safalt_fb200_deep40km_1Hz_p5`
# (133,699,789 tets) is the same mesh with those closed to ZERO, so swapping the
# parent fixes three quarters of the problem for free.
#
# It cannot reuse the shipped collar: measured by wall_compare.py, the 1 Hz
# parent's vertical wall is 46,041 triangles / 25,069 verts against the
# refine2/intermediate wall's 26,704 / 15,205 (its LEB pass bisected the wall).
# So the collar is rebuilt against the refined wall, same uniform 2,500 m recipe
# that beat the gate-driven collar 29x -- this is the BASE, which the LEB stage
# then refines to the gate.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
BC=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/.claude/worktrees/shakeoutbox-gate-close/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_intermediate/code
TMP=$A/meshing_shakeoutbox_gate/build_tmp
cd $TMP
$PY $BC/build_collar.py \
  --parent $A/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5 \
  --out $TMP/collar_1hz_uniform_2500.npz \
  --h-mode const --h-const 2500 --h-max 2500
