#!/bin/zsh
# Authoritative final numbers on the two shipped products.
CODE=${0:A:h}
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP

echo "########## INTERMEDIATE"
$PY $CODE/census_box.py --muscal --mesh $R/int_p4.puml.h5 --parent-tets 0 \
    --gate 0.6667 --gate 0.8 --dump $TMP/fail_int_SHIP.npz \
    > census_int_SHIP.log 2>&1
echo "  rc=$?"
$PY $CODE/residual_impact.py --dump $TMP/fail_int_SHIP.npz --gate 0.6667 \
    --order 3 --label "INTERMEDIATE ship, 0.5 Hz @ p3" > residual_int_SHIP.log 2>&1
echo "  rc=$?"

echo "########## HEAVY"
$PY $CODE/census_box.py --muscal --mesh $R/heavy_p3.puml.h5 --parent-tets 0 \
    --gate 0.6667 --gate 0.8 --dump $TMP/fail_heavy_SHIP.npz \
    > census_heavy_SHIP.log 2>&1
echo "  rc=$?"
$PY $CODE/residual_impact.py --dump $TMP/fail_heavy_SHIP.npz --gate 0.8 \
    --order 5 --label "HEAVY ship, 1 Hz @ p5" > residual_heavy_SHIP.log 2>&1
echo "  rc=$?"

echo "########## fault identity, both, vs the ORIGINAL shipped meshes"
$PY $CODE/check_fault_identity.py \
  --parent $A/meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5 \
  --new $R/int_p4.puml.h5 > fid_int_SHIP.log 2>&1
echo "  int rc=$?"
$PY $CODE/check_fault_identity.py \
  --parent $A/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5 \
  --new $R/heavy_p3.puml.h5 > fid_heavy_SHIP.log 2>&1
echo "  heavy rc=$?"
echo ALL DONE
