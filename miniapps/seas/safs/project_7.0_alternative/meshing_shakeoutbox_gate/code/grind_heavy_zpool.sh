#!/bin/zsh
# HEAVY -> zero cells below the 1 Hz @ p5 gate (Vs/dx >= 0.8), locally.
#
# THE COMBINATION THAT WAS NEVER TRIED HERE: --seed-on-gate WITH --pool zpool.
#   * --pool gate diverged on this mesh (369 -> 437, worst 0.494 -> 0.289): at
#     0.8 a bisected surface cell's barycentre reaches slow material faster than
#     dx shrinks.  The z-pooled lower bound is the cure for exactly that.
#   * mesh-wide zpool is unaffordable (5,287,420 cells flagged), but SEEDED on
#     the 369 measured failures the patch stays small -- the same combination
#     took the intermediate 124 -> 21.
#
# hops 110 FAILED THE OTHER WAY: a patch wide enough to stop rim-stalling also
# contains the pooled bound's GLOBAL demand -- 1.04M cells marked, +500k
# cells/ROUND. Fixed by --pool-radius 4: always refine measured gate failures,
# but apply the pooled bound only within 4 hops of one, which is where the
# treadmill actually bites.
#
# PASS 1 (hops 30) gave 369 -> 340 but stopped with "all 222 terminal edges
# frozen (rim 222)": the LEPP chains that would finish the job wander further
# than a 30-hop patch reaches. The failures occupy one 47x45 km corner, so a much
# wider patch is still cheap -- widened to 110 hops.
#
# GREEDY WITH ROLLBACK. Every pass is accepted only if the gate count STRICTLY
# decreases; otherwise the output is deleted and the previous mesh stands. This
# is deliberate: an earlier loop deleted its input unconditionally and threw away
# a better mesh (int_z1, 21 failures) in favour of a worse one (int_z2, 37).
set -e
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
CVM=/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc
TMP=$A/meshing_shakeoutbox_gate/build_tmp
R=$A/meshing_shakeoutbox_gate/results
cd $TMP

BEST=$R/hzp1.puml.h5
ORIG=$R/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5
BESTN=340
echo "start: $BESTN failing @0.8"

for i in 2 3 4 5 6 7 8; do
  OUT=$R/hzp$i.puml.h5
  # allow fault splitting from pass 3 on, once the volume-only route has had a
  # fair run; it is deck-visible (DR facet count) so it is not the first resort.
  EXTRA=""
  [ $i -ge 4 ] && EXTRA="--allow-fault-split"
  echo "=========== pass $i  (zpool r4, seeded, hops 60) $EXTRA"
  $PY $CODE/leb_muscal.py --mesh $BEST --cvm $CVM --out $OUT \
      --gate 0.8 --vs-source muscal --seed-on-gate --pool zpool $EXTRA \
      --hops 60 --pool-radius 4 --max-rounds 400 --stats $TMP/hzp${i}_stats.json \
      > close_heavy_zp$i.log 2>&1 || { echo "  leb FAILED"; break; }
  N=$($PY $CODE/census_box.py --muscal --mesh $OUT --parent-tets 0 --gate 0.8 \
        2>/dev/null | awk '/^    TOTAL/{gsub(",","",$2); print $2; exit}')
  T=$($PY -c "import h5py;print(h5py.File('$OUT','r')['connect'].shape[0])")
  echo "  pass $i -> $N failing, $T tets   (best so far $BESTN)"
  if [ -z "$N" ]; then echo "  census FAILED; keeping $BEST"; rm -f $OUT; break; fi
  if [ "$N" -ge "$BESTN" ]; then
    echo "  REJECTED (no strict improvement); keeping $BEST"
    rm -f $OUT
    break
  fi
  [ "$BEST" != "$ORIG" ] && rm -f $BEST      # never delete the shipped original
  BEST=$OUT; BESTN=$N
  [ "$N" = "0" ] && { echo "GATE CLOSED at pass $i"; break; }
done
echo "BEST $BEST  ($BESTN failing @0.8)"
