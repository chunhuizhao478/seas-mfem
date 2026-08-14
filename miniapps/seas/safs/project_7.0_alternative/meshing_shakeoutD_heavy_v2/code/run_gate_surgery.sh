#!/bin/zsh
# run_gate_surgery.sh -- batch gate surgery + finalize, unattended.
#   extract_multi (one pass) -> parallel refills with per-tile edge targets ->
#   stitch_multi (one rewrite) -> alt_finalize (census, LEB sweep, ASSERT 0,
#   checks, promote)
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
LOG=$B/gate_surgery.log
echo "=== gate surgery batch start $(date) ===" >> $LOG

$PY -u code/extract_multi.py --mesh $B/gate_state.mesh --boxes $B/gate_boxes.json \
    --outdir $B/tiles >> $LOG 2>&1 || { echo ABORT-extract; exit 1; }

# parallel refills, 4 at a time; each tile's edge target comes from its meta
ls $B/tiles/tile_*[0-9].msh | sed 's/\.msh$//' | xargs -P 4 -I{} zsh -c '
  PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
  base={}
  et=$($PY -c "import numpy as np;print(float(np.load(\"${base}_meta.npz\")[\"et\"]))")
  for seed in 0 1 2; do
    if $PY -u code/refill_box.py --box $base.msh --out ${base}_refill.msh \
         --shuffle $seed --edge-target $et --polish 2 >> $base.log 2>&1; then
      exit 0
    fi
  done
  echo "REFILL-FAILED $base" >> build_tmp/gate_surgery.log
'
nref=$(ls $B/tiles/tile_*_refill.msh 2>/dev/null | wc -l | tr -d ' ')
echo "refills done: $nref" >> $LOG

$PY -u code/stitch_multi.py --global-mesh $B/gate_state.mesh --tiledir $B/tiles \
    --out $B/v2_final4.mesh >> $LOG 2>&1 || { echo ABORT-stitch; exit 1; }

zsh code/alt_finalize.sh $B/v2_final4.mesh >> $LOG 2>&1
tail -1 $B/alt_finalize.log | tee -a $LOG
