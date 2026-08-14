#!/bin/zsh
# surgery_driver.sh -- clear the fault-locked sliver residue by cavity refill,
# box by box, chaining the global mesh state.
#
# Per box: extract (frozen 999 skin + fault, records removed tet ids) -> refill
# (tetgen -pq1.4/10Y on the box's own facets; retries --shuffle 1,2 on facet
# loss) -> stitch (asserts frozen-skin exact match + volume equality to 1e-9).
# Any assertion failure aborts the whole chain -- a bad splice must never
# propagate into the next box's extraction.
#
# Why refill and not splice-safe surgery: measured, 9 whole-mesh -optim passes
# flatlined at 435 and remove_sliver_tets found "NO splice-safe remesh" on every
# worst cluster -- but the refill of the first hotspot box produced eta_min
# 0.0972, ZERO below 0.05, with 0 of 5,588 facets lost. The residue is a local
# minimum of the old connectivity; replacing the pocket's interior escapes it.
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
CUR=${1:-$B/v2_s1.mesh}
LOG=$B/surgery.log

echo "=== surgery driver: start $(date)  base $CUR ===" >> $LOG

# 1. re-localize on the CURRENT mesh (earlier boxes may already be fixed)
$PY -u code/sliver_localize_medit.py $CUR >> $LOG 2>&1 || { echo ABORT-localize; exit 1; }
$PY -u code/make_surgery_boxes.py >> $LOG 2>&1 || { echo ABORT-boxes; exit 1; }

NB=$($PY -c "import json;print(len(json.load(open('$B/surgery_boxes.json'))))")
echo "boxes: $NB" >> $LOG

k=0
while [ $k -lt $NB ]; do
  read -r CX CY HX HY Z0 Z1 NS <<< $($PY -c "
import json; b=json.load(open('$B/surgery_boxes.json'))[$k]
print(b[0],b[1],b[2],b[3],b[4],b[5],b[6])")
  echo "--- box $k: $NS slivers at ($CX,$CY) ---" >> $LOG
  $PY -u code/extract_patch_medit.py --mesh $CUR --out $B/sb$k.msh \
      --meta $B/sb${k}_meta.npz --cx $CX --cy $CY --hx $HX --hy $HY \
      --z0 $Z0 --z1 $Z1 >> $LOG 2>&1 || { echo "ABORT-extract-$k"; exit 1; }
  ok=0
  for seed in 0 1 2; do
    if $PY -u code/refill_box.py --box $B/sb$k.msh --out $B/sb${k}_refill.msh \
         --shuffle $seed >> $LOG 2>&1; then ok=1; break; fi
    echo "  refill seed $seed failed, retrying" >> $LOG
  done
  [ $ok -eq 1 ] || { echo "ABORT-refill-$k"; exit 1; }
  NXT=$B/v2_srg$k.mesh
  $PY -u code/stitch_patch.py --global-mesh $CUR --patch $B/sb${k}_refill.msh \
      --meta $B/sb${k}_meta.npz --out $NXT >> $LOG 2>&1 || { echo "ABORT-stitch-$k"; exit 1; }
  # drop the previous intermediate (keep the original base)
  [ "$CUR" != "$B/v2_s1.mesh" ] && rm -f $CUR ${CUR%.mesh}.sol
  rm -f $B/sb$k.msh $B/sb${k}_refill.msh
  CUR=$NXT
  k=$((k+1))
done

echo "=== all boxes done: $CUR ===" >> $LOG
# final verification
$PY -u code/sliver_localize_medit.py $CUR >> $LOG 2>&1
mv $CUR $B/v2_final.mesh
echo "=== surgery complete $(date): $B/v2_final.mesh ===" >> $LOG
echo "SURGERY-DONE $B/v2_final.mesh"
