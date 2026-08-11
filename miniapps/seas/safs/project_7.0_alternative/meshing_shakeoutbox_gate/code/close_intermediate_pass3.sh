#!/bin/zsh
# Pass 3 -- the tail. 886 cells left, 99.7 % of them free-surface-resting, median
# refinement needed 1.10x and an ideal bill of 1,974 cells, so this is the
# convergence tail and not a structural floor.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP
/usr/bin/time -l $PY $CODE/leb_muscal.py \
  --mesh $R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal_p2.puml.h5 \
  --cvm $CVM \
  --out $R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal_final.puml.h5 \
  --gate 0.6667 --vs-source muscal --pool gate --hops 14 --max-rounds 90 \
  --stats $TMP/close_intermediate_p3_stats.json
