#!/bin/zsh
# INTERMEDIATE: close the 0.5 Hz @ p3 gate (Vs/dx >= 0.6667) on the MERGED mesh,
# judged against MUSCAL.
#
# Refining after the weld rather than the collar alone, because BOTH blocks have
# a residual on MUSCAL (parent 1,160 / collar 38,121) and because a merged-mesh
# pass has no frozen wall to stall against.  Fault edges stay frozen, so every
# fault-referenced deck input transfers unchanged.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
IN=$A/meshing_shakeoutbox_intermediate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox.puml.h5
OUT=$A/meshing_shakeoutbox_gate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5
cd $TMP
/usr/bin/time -l $PY $CODE/leb_muscal.py \
  --mesh $IN --cvm $CVM --out $OUT \
  --gate 0.6667 --vs-source muscal --pool gate --hops 4 --max-rounds 40 \
  --stats $TMP/close_intermediate_stats.json
