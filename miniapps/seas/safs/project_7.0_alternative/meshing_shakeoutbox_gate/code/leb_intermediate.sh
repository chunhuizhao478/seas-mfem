#!/bin/zsh
# Close the 0.5 Hz @ p3 gate (Vs/dx >= 0.6667) inside the INTERMEDIATE collar.
#
# Base = the shipped uniform 2,500 m collar, which is on the gate's COARSE
# branch: big cells put their barycentres in competent rock 800 m down.  That
# works for 99.7 % of them; the 42,228 that fail are the columns where even the
# coarse branch does not exist (basins), and no cell size between ~350 m and
# ~2,000 m passes there.  Bisection walks down through that forbidden zone, so
# the target must be the DOWNWARD-CLOSED one -- min Vs over the cell's own
# bounding box, which is a lower bound on every descendant and therefore makes
# compliance monotone.  --target bary or z-only pooling both treadmill here
# (measured: 15 -> 82 failures for +186 % tets).
#
# --all: the collar is standalone, so there is no outer mesh to protect and rim
# freezing can be dropped entirely.  Wall edges stay frozen regardless, so the
# weld interface is preserved exactly.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
BC=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/.claude/worktrees/shakeoutbox-gate-close/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_intermediate/code
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
cd $TMP
/usr/bin/time -l $PY $BC/leb_collar.py \
  --collar $A/meshing_shakeoutbox_intermediate/build_tmp/collar_uniform_2500.npz \
  --parent $A/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5 \
  --cvm $CVM \
  --out $TMP/collar_int_gate06667.npz \
  --gate 0.6667 --target bbox --all --max-rounds 40 \
  --stats $TMP/leb_intermediate_stats.json
