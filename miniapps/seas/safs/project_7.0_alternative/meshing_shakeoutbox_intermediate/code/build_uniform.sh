#!/bin/zsh
# Uniform far-field collar: no gate-driven surface refinement, no interior seed.
# The collar grades from the parent's own wall spacing (~309 m at the top rim)
# out to a uniform h, so the seam is a graded transition rather than a jump.
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
BASE=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
H=${1:-2500}
cd $CODE
$PY build_collar.py \
  --parent $BASE/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5 \
  --out $BASE/meshing_shakeoutbox_intermediate/build_tmp/collar_uniform_${H}.npz \
  --h-mode const --h-const $H --h-max $H
