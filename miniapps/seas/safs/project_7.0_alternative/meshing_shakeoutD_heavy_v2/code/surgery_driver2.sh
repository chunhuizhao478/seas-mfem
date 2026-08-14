#!/bin/zsh
# surgery_driver2.sh -- round 2 of cavity-refill surgery, threshold eta < 0.1.
#
# A separate file from surgery_driver.sh ON PURPOSE: round 1 was still executing
# when round 2 was prepared, and zsh reads scripts lazily -- editing a running
# script corrupts its execution mid-flight.
#
# Round 1 boxes only the sub-0.05 clusters (the localizer's old hardcoded
# threshold), so the 0.05-0.1 band survives it and pins eta_min at ~0.05. This
# round re-localizes at --thr 0.1 and refills those clusters too. Refilled
# pockets measure local eta_min ~0.097 against the fault's exact geometric
# ceiling of 0.106, so the projected outcome BEATS the small-domain benchmark
# (eta_min 0.0763, 127 below 0.1).
#
#   usage: surgery_driver2.sh <base.mesh> [thr] [final_name]
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
CUR=${1:?base mesh required}
THR=${2:-0.1}
FINAL=${3:-$B/v2_final2.mesh}
LOG=$B/surgery2.log
BASE0=$CUR

echo "=== surgery round 2: start $(date)  base $CUR  thr $THR ===" >> $LOG
$PY -u code/sliver_localize_medit.py $CUR --thr $THR >> $LOG 2>&1 || { echo ABORT-localize; exit 1; }
$PY -u code/make_surgery_boxes.py >> $LOG 2>&1 || { echo ABORT-boxes; exit 1; }
NB=$($PY -c "import json;print(len(json.load(open('$B/surgery_boxes.json'))))")
echo "boxes: $NB" >> $LOG

k=0
while [ $k -lt $NB ]; do
  read -r CX CY HX HY Z0 Z1 NS <<< $($PY -c "
import json; b=json.load(open('$B/surgery_boxes.json'))[$k]
print(b[0],b[1],b[2],b[3],b[4],b[5],b[6])")
  echo "--- r2 box $k: $NS slivers at ($CX,$CY) ---" >> $LOG
  # retry ladder: shuffle seeds, then GROW the box (+400 m per step -- a lost
  # skin facet usually lands interior to a larger box; a lost fault facet gets a
  # different local constellation). A box that still refuses is SKIPPED and
  # recorded, not allowed to kill the whole chain: its slivers stay, the census
  # at the end reports them honestly.
  ok=0
  for grow in 0 400 800; do
    $PY -u code/extract_patch_medit.py --mesh $CUR --out $B/r2b$k.msh \
        --meta $B/r2b${k}_meta.npz --cx $CX --cy $CY \
        --hx $((HX + grow)) --hy $((HY + grow)) \
        --z0 $((Z0 - grow)) --z1 $((Z1 + grow)) >> $LOG 2>&1 || { echo "ABORT-extract-$k"; exit 1; }
    for seed in 0 1 2 3; do
      if $PY -u code/refill_box.py --box $B/r2b$k.msh --out $B/r2b${k}_refill.msh \
           --shuffle $seed >> $LOG 2>&1; then ok=1; break; fi
      echo "  refill grow=$grow seed=$seed failed" >> $LOG
    done
    [ $ok -eq 1 ] && break
  done
  if [ $ok -ne 1 ]; then
    echo "SKIP-box-$k (refill unrecoverable after 3 grows x 4 seeds; slivers remain)" >> $LOG
    rm -f $B/r2b$k.msh $B/r2b${k}_refill.msh
    k=$((k+1))
    continue
  fi
  NXT=$B/v2r2_$k.mesh
  $PY -u code/stitch_patch.py --global-mesh $CUR --patch $B/r2b${k}_refill.msh \
      --meta $B/r2b${k}_meta.npz --out $NXT >> $LOG 2>&1 || { echo "ABORT-stitch-$k"; exit 1; }
  [ "$CUR" != "$BASE0" ] && rm -f $CUR ${CUR%.mesh}.sol
  rm -f $B/r2b$k.msh $B/r2b${k}_refill.msh
  CUR=$NXT
  k=$((k+1))
done

echo "=== r2 boxes done: $CUR ===" >> $LOG
$PY -u code/sliver_localize_medit.py $CUR --thr $THR >> $LOG 2>&1
mv $CUR $FINAL
echo "=== surgery round 2 complete $(date): $FINAL ===" >> $LOG
echo "SURGERY2-DONE $FINAL"
