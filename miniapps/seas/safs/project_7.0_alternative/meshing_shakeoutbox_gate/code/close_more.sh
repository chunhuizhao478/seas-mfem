#!/bin/zsh
# Keep grinding the tail on BOTH products. Each pass reseeds a fresh, wide patch
# on whatever is left; the count falls ~4x per pass for ~10x the rounds, and the
# cells are cheap (median refinement still needed 1.10-1.20x).
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP

echo "########## INTERMEDIATE pass 4"
$PY $CODE/leb_muscal.py \
  --mesh $R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal_final.puml.h5 \
  --cvm $CVM --out $R/int_p4.puml.h5 \
  --gate 0.6667 --vs-source muscal --pool gate --hops 20 --max-rounds 140 \
  --stats $TMP/int_p4_stats.json > close_int_p4.log 2>&1
echo "  rc=$?"

echo "########## HEAVY pass 3"
$PY $CODE/leb_muscal.py \
  --mesh $R/safalt_fb200_1Hz_p5_shakeoutbox_muscal_final.puml.h5 \
  --cvm $CVM --out $R/heavy_p3.puml.h5 \
  --gate 0.8 --vs-source muscal --pool gate --hops 20 --max-rounds 140 \
  --stats $TMP/heavy_p3_stats.json > close_heavy_p3.log 2>&1
echo "  rc=$?"
echo ALL DONE
