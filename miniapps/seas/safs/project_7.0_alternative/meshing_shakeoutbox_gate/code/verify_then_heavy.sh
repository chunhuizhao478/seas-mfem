#!/bin/zsh
# 1) verify the finished intermediate, 2) start the heavy by welding the new
# collar onto the 1 Hz @ p5 parent, 3) census the weld.
CODE=${0:A:h}
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
INT=$R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal_p2.puml.h5
cd $TMP

echo "########## 1. INTERMEDIATE census on MUSCAL"
$PY $CODE/census_box.py --muscal --mesh $INT --parent-tets 0 \
    --gate 0.6667 --gate 0.8 --dump $TMP/fail_int_final_muscal.npz \
    > census_int_final_muscal.log 2>&1
echo "  rc=$?"

echo "########## 2. INTERMEDIATE fault identity vs the shipped mesh"
$PY $CODE/check_fault_identity.py \
    --parent $A/meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5 \
    --new $INT > fault_identity_int.log 2>&1
echo "  rc=$?"

echo "########## 3. HEAVY merge (parent swap to the 1 Hz p5 mesh)"
$CODE/merge_heavy.sh > merge_heavy.log 2>&1
echo "  rc=$?"

echo "########## 4. HEAVY census on MUSCAL"
$PY $CODE/census_box.py --muscal \
    --mesh $TMP/safalt_fb200_1Hz_p5_shakeoutbox_merged.puml.h5 \
    --parent-tets 133699789 --gate 0.6667 --gate 0.8 \
    --dump $TMP/fail_heavy_merged_muscal.npz \
    > census_heavy_merged_muscal.log 2>&1
echo "  rc=$?"
echo ALL DONE
