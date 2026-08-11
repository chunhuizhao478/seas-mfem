#!/bin/zsh
# The true baseline: every mesh judged against MUSCAL, the SOURCE velocity
# model, rather than the deck's 250 m-binned resample.
CODE=${0:A:h}
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
TMP=$A/meshing_shakeoutbox_gate/build_tmp
cd $TMP

echo "########## INTERMEDIATE product (parent 15,979,903 + collar)"
$PY $CODE/census_box.py --muscal \
  --mesh $A/meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5 \
  --parent-tets 15979903 --gate 0.6667 --gate 0.8 \
  --dump $TMP/fail_intermediate_muscal.npz > census_intermediate_muscal.log 2>&1
echo "rc=$?"

echo "########## 1 Hz p5 parent (the heavy mesh's parent)"
$PY $CODE/census_box.py --muscal \
  --mesh $A/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5 \
  --parent-tets 0 --gate 0.6667 --gate 0.8 \
  --dump $TMP/fail_parent1hz_muscal.npz > census_parent1hz_muscal.log 2>&1
echo "rc=$?"
echo ALL DONE
