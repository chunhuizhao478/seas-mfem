#!/bin/zsh
# Tail pass: seed on the MEASURED failures (tiny patch), drive with the POOLED
# bound (kills the oscillation the gate target cannot).
#   $1 = int | heavy
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP
case $1 in
  int)   IN=$R/int_p4.puml.h5;   OUT=$R/int_p5.puml.h5;   G=0.6667 ;;
  heavy) IN=$R/heavy_p3.puml.h5; OUT=$R/heavy_p4.puml.h5; G=0.8 ;;
  *) echo "usage: $0 {int|heavy}"; exit 2 ;;
esac
/usr/bin/time -l $PY $CODE/leb_muscal.py \
  --mesh $IN --cvm $CVM --out $OUT \
  --gate $G --vs-source muscal --seed-on-gate --pool half \
  --hops 8 --max-rounds 60 --stats $TMP/tail_$1_stats.json
