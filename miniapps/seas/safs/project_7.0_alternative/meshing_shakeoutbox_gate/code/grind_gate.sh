#!/bin/zsh
# Iterated GATE-target passes, to zero.
#
# Why this beats the pooled bound for the last cells: the pooled bound is the
# convergence GUARANTEE, but it wants refinement everywhere a free-surface cell
# is coarser than Vs(depth 0)/gate -- 3,182,350 cells mesh-wide, a 78 % patch.
# Bounding that with a seeded patch just moves the problem to the patch RIM
# (measured: grind 1 stalled with 1 terminal edge and 3,898 frozen at the rim).
#
# The gate target has no such appetite: it touches only cells that actually
# fail. It cannot oscillate forever either -- the size a cell must reach is
# bounded below by Vs(depth 0)/gate > 0, because MUSCAL has nothing shallower
# than the surface -- so each pass strictly reduces the count and the process
# terminates. Measured: pass 4 took 246 -> 124 for +45,491 tets. Iterate it.
#
#   $1 = starting mesh   $2 = tag
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP

CUR=$1
TAG=${2:-g}
PREV=999999999
for i in 1 2 3 4 5 6 7 8 9 10 11 12; do
  OUT=$R/${TAG}$i.puml.h5
  echo "=========== $TAG$i : $(basename $CUR) -> $(basename $OUT)"
  $PY $CODE/leb_muscal.py --mesh $CUR --cvm $CVM --out $OUT \
      --gate 0.6667 --vs-source muscal --seed-on-gate --pool gate \
      --hops 30 --max-rounds 150 --stats $TMP/${TAG}${i}_stats.json \
      > close_${TAG}$i.log 2>&1 || { echo "  leb FAILED"; break; }
  N=$($PY $CODE/census_box.py --muscal --mesh $OUT --parent-tets 0 --gate 0.6667 \
        2>/dev/null | awk '/^    TOTAL/{gsub(",","",$2); print $2; exit}')
  T=$($PY -c "import h5py;print(h5py.File('$OUT','r')['connect'].shape[0])")
  echo "  after $TAG$i : $N failing, $T tets"
  [ "$CUR" != "$1" ] && rm -f $CUR
  CUR=$OUT
  if [ "$N" = "0" ]; then echo "GATE CLOSED at $TAG$i"; break; fi
  if [ "$N" -ge "$PREV" ]; then echo "NO PROGRESS ($PREV -> $N); stopping"; break; fi
  PREV=$N
done
echo "FINAL $CUR"
