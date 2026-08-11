#!/bin/zsh
# HEAVY step 1: weld the new uniform collar onto the 1 Hz @ p5 parent.
#
# This is the parent SWAP.  The shipped heavy mesh welds its collar onto
# `safalt_fb200_deep40km_refine2` (122,162,105 tets), which was only ever
# certified at 0.5 Hz/p3 and fails the 1 Hz gate on 274,299 of its own cells --
# 75 % of that product's total.  `safalt_fb200_deep40km_1Hz_p5` (133,699,789)
# is the same mesh with those closed, and measures 0 failures at BOTH gates on
# the deck nc and 7,648 at 0.8 on MUSCAL.
#
# The collar had to be rebuilt because the 1 Hz parent's vertical wall is
# 46,041 triangles / 25,069 verts against refine2's 26,704 / 15,205 -- its LEB
# pass bisected the wall, so the two are NOT interchangeable (wall_compare.py).
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
BC=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/.claude/worktrees/shakeoutbox-gate-close/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_intermediate/code
TMP=$A/meshing_shakeoutbox_gate/build_tmp
P=$A/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5
OUT=$TMP/safalt_fb200_1Hz_p5_shakeoutbox_merged.puml.h5
cd $TMP
/usr/bin/time -l $PY $BC/merge_collar.py \
  --parent $P --collar $TMP/collar_1hz_uniform_2500.npz --out $OUT
