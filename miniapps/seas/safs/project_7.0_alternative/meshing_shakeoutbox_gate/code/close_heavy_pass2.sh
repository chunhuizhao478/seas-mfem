#!/bin/zsh
# HEAVY pass 2. Pass 1 ended at 7,462 marked with 390 terminal edges frozen at
# its 5-hop patch RIM; reseeding a fresh, wider patch on the residual is what
# releases them.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP
/usr/bin/time -l $PY $CODE/leb_muscal.py \
  --mesh $R/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5 \
  --cvm $CVM \
  --out $R/safalt_fb200_1Hz_p5_shakeoutbox_muscal_final.puml.h5 \
  --gate 0.8 --vs-source muscal --pool gate --hops 12 --max-rounds 80 \
  --stats $TMP/close_heavy_p2_stats.json
