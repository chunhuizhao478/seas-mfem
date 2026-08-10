#!/bin/zsh
# Weld the uniform far-field collar onto BOTH deep ALT parents (shared wall).
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
WT=${0:A:h}/..
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
COL=$A/meshing_shakeoutbox_intermediate/build_tmp/collar_uniform_2500.npz

echo "=========== INTERMEDIATE ==========="
cd ${0:A:h}
$PY merge_collar.py \
  --parent $A/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5 \
  --collar $COL \
  --out $A/meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5

echo "=========== HEAVY ==========="
cd $A/../../../../.claude/worktrees/shakeoutbox-alt-extend/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_heavy/code 2>/dev/null || cd ${0:A:h}
$PY merge_collar.py \
  --parent $A/meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5 \
  --collar $COL \
  --out $A/meshing_shakeoutbox_heavy/results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5
