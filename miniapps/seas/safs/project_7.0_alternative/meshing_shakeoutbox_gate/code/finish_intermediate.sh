#!/bin/zsh
# Weld the gate-closed collar onto the intermediate parent, then verify.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CODE=${0:A:h}
BC=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/.claude/worktrees/shakeoutbox-gate-close/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_intermediate/code
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
OUT=$A/meshing_shakeoutbox_gate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox_gate.puml.h5
P=$A/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5
cd $TMP

echo "########## MERGE"
$PY $BC/merge_collar.py --parent $P --collar $TMP/collar_int_gate06667.npz --out $OUT

echo "########## VERIFY (acceptance gates P1-P9)"
$PY $BC/verify_extension_lean.py --parent $P --new $OUT --cvm $CVM --gate 0.6667 || true

echo "########## CENSUS"
NP=$($PY -c "import h5py,sys; print(h5py.File('$P','r')['connect'].shape[0])")
$PY $CODE/census_box.py --mesh $OUT --parent-tets $NP --gate 0.6667 --gate 0.8 \
    --dump $TMP/fail_intermediate_gate.npz
