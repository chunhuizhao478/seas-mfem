#!/bin/zsh
# Drive the INTERMEDIATE mesh to ZERO cells below 0.5 Hz at p3.
#
# The gate-as-target floors out around 100 cells and oscillates, because
# bisecting a surface cell moves its barycentre nearer z=0 where Vs is slower.
# --pool zpool is the exact cure: min Vs over the cell's own VERTICAL extent is
# a true lower bound on every descendant (Vs bottoms out at depth 0), so
# compliance becomes monotone. Unlike bbox/half it costs 5 lookups and no
# minimum_filter, so it is affordable mesh-wide.
#
#   $1 = input mesh, $2 = output mesh, $3 = hops, $4 = max-rounds
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP
/usr/bin/time -l $PY $CODE/leb_muscal.py \
  --mesh ${1} --cvm $CVM --out ${2} \
  --gate 0.6667 --vs-source muscal --pool zpool \
  --hops ${3:-8} --max-rounds ${4:-120} \
  --stats $TMP/$(basename ${2} .puml.h5)_stats.json
