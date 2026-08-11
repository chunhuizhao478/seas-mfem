#!/bin/zsh
# HEAVY step 2: close the 1 Hz @ p5 gate (Vs/dx >= 0.8) on the merged mesh.
#
# 0.8 is the binding gate -- the heavy mesh must resolve 0.5 Hz at p3 (0.6667)
# AND 1 Hz at p5 (0.8000), and the second implies the first.
#
# --pool gate: on MUSCAL the refinement target can BE the gate.  The pooled
# lower bound exists to stop a treadmill, and the treadmill is a property of the
# deck nc's 250 m binning (Vs steps 3.78x across the z=-125 m bin edge, so one
# bisection could drop a child's Vs 4x).  MUSCAL's ladder is 50 m through the
# top 500 m and Vs BOTTOMS OUT at the surface rather than falling off a cliff,
# so refining converges on its own -- measured on the intermediate mesh, where
# the pooled target flagged 587,899-1,456,239 cells against 39,281 that
# actually fail, and the gate target closed 39,281 -> 886 for +17.9 % tets.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP
/usr/bin/time -l $PY $CODE/leb_muscal.py \
  --mesh $TMP/safalt_fb200_1Hz_p5_shakeoutbox_merged.puml.h5 \
  --cvm $CVM \
  --out $R/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5 \
  --gate 0.8 --vs-source muscal --pool gate --hops 5 --max-rounds 55 \
  --stats $TMP/close_heavy_stats.json
