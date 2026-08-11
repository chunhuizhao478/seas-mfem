#!/bin/zsh
# Drive the INTERMEDIATE mesh to ZERO cells below 0.5 Hz at p3.
#
# Strategy: seed the patch on the MEASURED failures only (a few hundred cells,
# so the patch stays ~1 % of the mesh) and drive it with the Z-POOLED lower
# bound, which makes compliance monotone inside that patch.  Mesh-wide zpool is
# the correct bound but flags 3,182,350 cells and a 78 % patch, because it asks
# every free-surface cell in the collar to be sized for Vs(depth 0) -- the
# expensive fine branch.  Seeding on failures buys the convergence guarantee
# only where the treadmill actually bites.
#
# Iterated: each round re-measures, so a pass that ends with survivors is simply
# followed by another.  Stops when a census reports 0, or when a pass fails to
# reduce the count (then the log is kept for diagnosis).
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP

CUR=$R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5
for i in 1 2 3 4 5 6 7 8; do
  OUT=$R/int_z$i.puml.h5
  echo "=========== grind $i : $(basename $CUR) -> $(basename $OUT)"
  $PY $CODE/leb_muscal.py --mesh $CUR --cvm $CVM --out $OUT \
      --gate 0.6667 --vs-source muscal --seed-on-gate --pool zpool \
      --hops 25 --max-rounds 200 --stats $TMP/int_z${i}_stats.json \
      > close_int_z$i.log 2>&1
  rc=$?
  echo "  leb rc=$rc"
  [ $rc -ne 0 ] && break
  N=$($PY $CODE/census_box.py --muscal --mesh $OUT --parent-tets 0 --gate 0.6667 \
        2>/dev/null | awk '/^    TOTAL/{gsub(",","",$2); print $2; exit}')
  T=$($PY -c "import h5py;print(h5py.File('$OUT','r')['connect'].shape[0])")
  echo "  after grind $i : $N failing, $T tets"
  [ -f "$CUR" ] && [ "$CUR" != "$R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5" ] && rm -f $CUR
  CUR=$OUT
  if [ "$N" = "0" ]; then echo "GATE CLOSED at grind $i"; break; fi
done
echo "FINAL $CUR"
