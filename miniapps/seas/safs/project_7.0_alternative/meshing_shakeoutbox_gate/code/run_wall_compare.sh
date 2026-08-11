#!/bin/zsh
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
cd $A/meshing_shakeoutbox_gate/build_tmp
$PY $CODE/wall_compare.py \
  --mesh $A/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5 \
  --mesh $A/meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5 \
  --mesh $A/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5
