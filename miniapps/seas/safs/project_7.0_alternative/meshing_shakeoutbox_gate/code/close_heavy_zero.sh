#!/bin/zsh
# Drive the HEAVY mesh to ZERO cells below the 1 Hz @ p5 gate (Vs/dx >= 0.8).
#
# Applies the two lessons that closed the intermediate, neither of which had been
# tried on the heavy:
#   A) --seed-on-gate --hops 30 : seed the patch on the MEASURED failures so it
#      stays a small % of the mesh even at 30 hops. The earlier heavy passes used
#      --hops 5/12 and stalled against their own patch RIM, which reads as a
#      convergence floor but is not one (intermediate: 37 -> 4 for +497 tets).
#   B) --allow-fault-split : release the fault-edge-pinned class, where a cell's
#      LONGEST edge is an edge of the fault's own triangulation and no volume
#      refinement can shorten it. Rivara inserts the edge MIDPOINT, which for a
#      fault edge lies exactly on the planar fault triangles sharing it, so the
#      fault surface and AREA are preserved bit-for-bit; only the triangulation
#      refines. On the intermediate this cost 63 tets and moved 9 of 160,280
#      facets. It IS deck-visible via the DR facet count.
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP

CUR=$R/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5
for i in 1 2 3 4 5; do
  OUT=$R/hz$i.puml.h5
  # fault splitting is allowed from round 2 on: round 1 shows whether the
  # residual is reachable without touching the fault at all.
  EXTRA=""
  [ $i -ge 2 ] && EXTRA="--allow-fault-split"
  echo "=========== heavy grind $i  $EXTRA"
  $PY $CODE/leb_muscal.py --mesh $CUR --cvm $CVM --out $OUT \
      --gate 0.8 --vs-source muscal --seed-on-gate --pool gate $EXTRA \
      --hops 30 --max-rounds 150 --stats $TMP/hz${i}_stats.json \
      > close_heavy_z$i.log 2>&1 || { echo "  leb FAILED"; break; }
  N=$($PY $CODE/census_box.py --muscal --mesh $OUT --parent-tets 0 --gate 0.8 \
        2>/dev/null | awk '/^    TOTAL/{gsub(",","",$2); print $2; exit}')
  T=$($PY -c "import h5py;print(h5py.File('$OUT','r')['connect'].shape[0])")
  echo "  after grind $i : $N failing @0.8, $T tets"
  [ "$CUR" != "$R/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5" ] && rm -f $CUR
  CUR=$OUT
  [ "$N" = "0" ] && { echo "GATE CLOSED at grind $i"; break; }
done
echo "FINAL $CUR"
