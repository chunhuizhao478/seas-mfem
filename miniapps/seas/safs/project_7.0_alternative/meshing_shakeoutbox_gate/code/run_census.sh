#!/bin/zsh
# Gate census on every mesh that matters, at both gates.
#   $1 = which  (intermediate | heavy | parent1hz | parentref2 | parentint)
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
TMP=$A/meshing_shakeoutbox_gate/build_tmp
cd $TMP

case $1 in
  intermediate)
    $PY $CODE/census_box.py \
      --mesh $A/meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5 \
      --parent-tets 15979903 --gate 0.6667 --gate 0.8 \
      --dump $TMP/fail_intermediate.npz ;;
  heavy)
    $PY $CODE/census_box.py \
      --mesh $A/meshing_shakeoutbox_heavy/results/safalt_fb200_deep40km_refine2_shakeoutbox.puml.h5 \
      --parent-tets 122162105 --gate 0.6667 --gate 0.8 \
      --dump $TMP/fail_heavy.npz ;;
  parentint)
    $PY $CODE/census_box.py \
      --mesh $A/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5 \
      --parent-tets 0 --gate 0.6667 --gate 0.8 ;;
  parentref2)
    $PY $CODE/census_box.py \
      --mesh $A/meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5 \
      --parent-tets 0 --gate 0.6667 --gate 0.8 ;;
  parent1hz)
    $PY $CODE/census_box.py \
      --mesh $A/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5 \
      --parent-tets 0 --gate 0.6667 --gate 0.8 ;;
  *) echo "usage: $0 {intermediate|heavy|parentint|parentref2|parent1hz}"; exit 2 ;;
esac
