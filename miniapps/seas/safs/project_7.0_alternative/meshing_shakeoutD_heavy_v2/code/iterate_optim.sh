#!/bin/zsh
# iterate_optim.sh -- the iterated-optim de-sliver campaign, as a driver.
#
# Encodes what was hand-driven on ALT (nine passes, 1,017,790 -> 435): per pass,
# run `mmg3d -opnbdy -optim -nosurf` with insertion allowed, measure the output,
# and continue while a pass removes >= STOPFRAC of the remaining sub-0.05 tets.
# Accept a pass only if the surface triangle count is unchanged (fault frozen).
#
#   usage: iterate_optim.sh <workdir> <input.mesh|input.puml.h5> <prefix> [max] [memMB]
#
# The same-basename .sol sidecar is deleted every pass -- mmg auto-loads it as an
# input metric and then refuses -optim (measured trap).
set -u
WD=$1; IN=$2; PFX=$3; MAXP=${4:-8}; MEM=${5:-24000}
STOPFRAC=0.30
CODE=${0:a:h}
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
MMG=/Users/chunhuizhao/miniforge/envs/mmg/bin/mmg3d_O3
cd $WD
LOG=${PFX}_campaign.log
touch ${PFX}.running

if [[ $IN == *.puml.h5 ]]; then
  echo "[convert] $IN -> ${PFX}_p0.mesh" >> $LOG
  $PY -u $CODE/puml_to_medit.py $IN ${PFX}_p0.mesh >> $LOG 2>&1 || { echo ABORT-convert; rm -f ${PFX}.running; exit 1; }
  CUR=${PFX}_p0.mesh
else
  CUR=$IN
fi

cd_line=$(cd $CODE/.. && pwd)
meas () { (cd $CODE/.. && $PY -u code/measure_eta_full.py $WD/$1) }
M0=$(meas $CUR | tee -a $LOG)
prev=$(echo $M0 | sed 's/.*lt005=\([0-9]*\).*/\1/')
tris0=$(echo $M0 | sed 's/.*tris=\([0-9]*\).*/\1/')
echo "[start] lt005=$prev tris=$tris0" >> $LOG

k=1
while [ $k -le $MAXP ]; do
  NXT=${PFX}_p$k.mesh
  rm -f ${CUR%.mesh}.sol
  echo "--- pass $k: $CUR -> $NXT $(date) ---" >> $LOG
  $MMG -in $CUR -out $NXT -opnbdy -optim -nosurf -hmin 20 -hmax 200000 \
       -hgrad 3 -m $MEM -v 2 >> $LOG 2>&1 || { echo "ABORT-mmg-$k" >> $LOG; rm -f ${PFX}.running; exit 1; }
  M=$(meas $NXT | tee -a $LOG)
  cur=$(echo $M | sed 's/.*lt005=\([0-9]*\).*/\1/')
  tris=$(echo $M | sed 's/.*tris=\([0-9]*\).*/\1/')
  if [ "$tris" != "$tris0" ]; then
    echo "ABORT: surface count changed $tris0 -> $tris (pass $k REJECTED, keeping $CUR)" >> $LOG
    rm -f $NXT ${NXT%.mesh}.sol ${PFX}.running
    exit 1
  fi
  # previous intermediate no longer needed (keep pass-0 conversion as the base)
  [ "$CUR" != "${PFX}_p0.mesh" ] && [ "$CUR" != "$IN" ] && rm -f $CUR ${CUR%.mesh}.sol
  CUR=$NXT
  removed=$((prev - cur))
  echo "[pass $k] lt005 $prev -> $cur (removed $removed)" >> $LOG
  # stop when the pass removed < STOPFRAC of what remained
  stop=$($PY -c "print(1 if $prev==0 or ($removed)/max($prev,1) < $STOPFRAC else 0)")
  prev=$cur
  [ "$stop" = "1" ] && break
  k=$((k+1))
done
mv $CUR ${PFX}_final.mesh
rm -f ${CUR%.mesh}.sol ${PFX}.running
echo "OPTIM-CAMPAIGN-DONE ${PFX}_final.mesh lt005=$prev" | tee -a $LOG
