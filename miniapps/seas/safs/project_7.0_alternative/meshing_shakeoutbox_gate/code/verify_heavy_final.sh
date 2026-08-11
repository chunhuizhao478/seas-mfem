#!/bin/zsh
# Final acceptance on the heavy product: the gate it must hold, and proof that
# the fault / free surface / BC contract survived the bisection.
CODE=${0:A:h}
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
H=$R/safalt_fb200_1Hz_p5_shakeoutbox_muscal_final.puml.h5
cd $TMP

echo "########## heavy census on MUSCAL, both gates"
$PY $CODE/census_box.py --muscal --mesh $H --parent-tets 0 \
    --gate 0.6667 --gate 0.8 --dump $TMP/fail_heavy_final_muscal.npz \
    > census_heavy_final_muscal.log 2>&1
echo "  rc=$?"

echo "########## heavy fault identity vs the merged (pre-refinement) mesh"
$PY $CODE/check_fault_identity.py \
    --parent $TMP/safalt_fb200_1Hz_p5_shakeoutbox_merged.puml.h5 \
    --new $H > fault_identity_heavy.log 2>&1
echo "  rc=$?"
echo ALL DONE
