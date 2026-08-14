#!/bin/zsh
# surgery_driver_gate.sh -- cavity-refill surgery on GATE failures.
#
# Same chain as the eta surgery (extract -> refill -> stitch, all assertions
# live) with one addition: each box carries an EDGE TARGET from gate_boxes.py
# (min Vs among its failures / 0.8 * 0.85), which refill_box turns into a tetgen
# -a volume cap and a polish -hmax. The fresh interior cannot recreate the gate
# failure -- unlike bisection, which fought the slow-basin CVM treadmill at ~45
# cells/round (measured) and stalled.
#
#   usage: surgery_driver_gate.sh <base.mesh> [final_name]
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
CUR=${1:?base mesh required}
FINAL=${2:-$B/v2_final4.mesh}
LOG=$B/surgery_gate.log
BASE0=$CUR

NB=$($PY -c "import json;print(len(json.load(open('$B/gate_boxes.json'))))")
echo "=== gate surgery: start $(date)  base $CUR  boxes $NB ===" >> $LOG

k=0
while [ $k -lt $NB ]; do
  read -r CX CY HX HY Z0 Z1 NS ET <<< $($PY -c "
import json; b=json.load(open('$B/gate_boxes.json'))[$k]
print(b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7])")
  echo "--- gate box $k: $NS fails at ($CX,$CY)  edge<=$ET ---" >> $LOG
  ok=0
  for grow in 0 400 800; do
    $PY -u code/extract_patch_medit.py --mesh $CUR --out $B/gb$k.msh \
        --meta $B/gb${k}_meta.npz --cx $CX --cy $CY \
        --hx $((HX + grow)) --hy $((HY + grow)) \
        --z0 $((Z0 - grow)) --z1 $Z1 >> $LOG 2>&1 || { echo "ABORT-extract-$k"; exit 1; }
    for seed in 0 1 2 3; do
      if $PY -u code/refill_box.py --box $B/gb$k.msh --out $B/gb${k}_refill.msh \
           --shuffle $seed --edge-target $ET --polish 2 >> $LOG 2>&1; then ok=1; break; fi
      echo "  refill grow=$grow seed=$seed failed" >> $LOG
    done
    [ $ok -eq 1 ] && break
  done
  if [ $ok -ne 1 ]; then
    echo "SKIP-box-$k (refill unrecoverable; gate fails remain here)" >> $LOG
    rm -f $B/gb$k.msh $B/gb${k}_refill.msh
    k=$((k+1))
    continue
  fi
  NXT=$B/vg_$k.mesh
  $PY -u code/stitch_patch.py --global-mesh $CUR --patch $B/gb${k}_refill.msh \
      --meta $B/gb${k}_meta.npz --out $NXT >> $LOG 2>&1 || { echo "ABORT-stitch-$k"; exit 1; }
  [ "$CUR" != "$BASE0" ] && rm -f $CUR ${CUR%.mesh}.sol
  rm -f $B/gb$k.msh $B/gb${k}_refill.msh
  CUR=$NXT
  k=$((k+1))
done

mv $CUR $FINAL
echo "=== gate surgery complete $(date): $FINAL ===" >> $LOG
echo "GATE-SURGERY-DONE $FINAL"
