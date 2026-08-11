#!/bin/zsh
# Close the 1 Hz @ p5 gate (Vs/dx >= 0.8) inside the HEAVY collar.
#
# Base = collar_1hz_uniform_2500.npz, built against the 1 Hz parent's REFINED
# wall (46,041 tris, vs 26,704 on refine2/intermediate -- the two are not
# interchangeable, measured by wall_compare.py).
#
# 0.8 is the binding gate for the heavy mesh: it must resolve 0.5 Hz at p3
# (0.6667) AND 1 Hz at p5 (0.8000), and the second implies the first.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
BC=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/.claude/worktrees/shakeoutbox-gate-close/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_intermediate/code
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
cd $TMP
/usr/bin/time -l $PY $BC/leb_collar.py \
  --collar $TMP/collar_1hz_uniform_2500.npz \
  --parent $A/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5 \
  --cvm $CVM \
  --out $TMP/collar_1hz_gate08.npz \
  --gate 0.8 --target bbox --all --max-rounds 40 \
  --stats $TMP/leb_heavy_stats.json
