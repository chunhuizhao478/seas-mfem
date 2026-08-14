#!/bin/zsh
# alt_finalize.sh -- from the surgery output to the DELIVERABLE:
#   the ALT ShakeOut-D heavy mesh resolving 1 Hz at p5 EVERYWHERE.
#
#   1. MEDIT -> PUML (medit_to_puml: derived offsets, interior-crack invariant)
#   2. gate census on NATIVE MUSCAL (never the deck cube -- 27x trap, measured)
#   3. LEB gate close to ZERO (max-rounds 250; only splits, measured to leave
#      eta bit-identical -- the quality won by surgery survives this stage)
#   4. re-census: assert 0 below gate
#   5. PUML structure checks (P1-P4) + promote to results/
#
# MEMORY POLICY, per the user's priority (ALT delivers first): the LEB stage
# needs ~18-20 GB and cannot share the box with a PREF mmg pass (~11-14 GB).
# Before stage 3 this script STOPS the PREF optim campaign cleanly (it resumes
# from its last accepted pass: iterate_optim.sh <wd> <last_pN.mesh> pref_optim).
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
IN=${1:-$B/v2_final2.mesh}
PUML0=$B/shakeoutD_alt_heavy_v2.puml.h5
GATED=$B/shakeoutD_alt_heavy_v2_gate.puml.h5
FINAL=results/safalt_shakeoutD_heavy_1Hz_p5_v2.puml.h5
CVM=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/seisol_quakeworx/safs_seisol_v3_2_0_RSSRW_ALTv2_freesurf/safs_material_cvm.nc
MUSCAL=/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc
LOG=$B/alt_finalize.log
PREF_MARK=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_preferred/meshing_shakeoutD_heavy/build_tmp/pref_optim.running

[ -f $IN ] || { echo "MISSING $IN -- surgery round 2 not finished?"; exit 1; }
mkdir -p results
echo "=== alt_finalize start $(date)  in: $IN ===" >> $LOG

echo "=== 1/5 MEDIT -> PUML ===" >> $LOG
$PY -u code/medit_to_puml.py $IN $PUML0 >> $LOG 2>&1 || { echo ABORT-convert; exit 1; }

echo "=== 2/5 gate census (native MUSCAL) BEFORE close ===" >> $LOG
$PY -u - $PUML0 $CVM $MUSCAL >> $LOG 2>&1 <<'PYEOF' || { echo ABORT-census; exit 1; }
import sys, numpy as np, h5py
sys.path.insert(0, 'code')
from material import Material
BOX=(132000.0,724000.0,3490000.0,4055000.0)
PAIRS=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
mat=Material(sys.argv[2],sys.argv[3],box=BOX,source='muscal',verbose=False)
f=h5py.File(sys.argv[1]); V=f['geometry'][:]; C=f['connect']; NT=C.shape[0]
nb=0; worst=9e9
for s in range(0,NT,6_000_000):
    t=C[s:s+6_000_000][:].astype(np.int64); p=V[t]
    dx=np.max(np.stack([np.linalg.norm(p[:,j]-p[:,i],axis=1) for i,j in PAIRS],1),1)
    vb=mat.at(p.mean(1)); r=np.where(vb>0,vb/dx,0.0)
    nb+=int((r<0.8).sum()); worst=min(worst,float(r.min())); del t,p,dx,vb,r
print(f'CENSUS {sys.argv[1]}: below-gate {nb:,} worst {worst:.4f}')
PYEOF

echo "=== 3/5 LEB gate close to ZERO (stopping PREF campaign for the slot) ===" >> $LOG
if [ -f $PREF_MARK ]; then
  echo "[slot] stopping PREF optim campaign (resumable from last accepted pass)" >> $LOG
  pkill -f iterate_optim.sh 2>/dev/null; sleep 1; pkill -x mmg3d_O3 2>/dev/null
  rm -f $PREF_MARK
fi
$PY -u code/leb_gate_close.py --mesh $PUML0 --parent $PUML0 --include-parent \
  --cvm $CVM --muscal $MUSCAL --vs-source muscal --gate 0.8 \
  --drive pooled --pool-frac 0.25 --max-rounds 250 --checkpoint-every 10 \
  --stats $B/alt_gate_stats.json --out $GATED >> $LOG 2>&1 || { echo ABORT-leb; exit 1; }

echo "=== 4/5 re-census: must be ZERO ===" >> $LOG
$PY -u - $GATED $CVM $MUSCAL >> $LOG 2>&1 <<'PYEOF' || { echo ABORT-census2; exit 1; }
import sys, numpy as np, h5py
sys.path.insert(0, 'code')
from material import Material
BOX=(132000.0,724000.0,3490000.0,4055000.0)
PAIRS=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
mat=Material(sys.argv[2],sys.argv[3],box=BOX,source='muscal',verbose=False)
f=h5py.File(sys.argv[1]); V=f['geometry'][:]; C=f['connect']; NT=C.shape[0]
nb=0; worst=9e9
for s in range(0,NT,6_000_000):
    t=C[s:s+6_000_000][:].astype(np.int64); p=V[t]
    dx=np.max(np.stack([np.linalg.norm(p[:,j]-p[:,i],axis=1) for i,j in PAIRS],1),1)
    vb=mat.at(p.mean(1)); r=np.where(vb>0,vb/dx,0.0)
    nb+=int((r<0.8).sum()); worst=min(worst,float(r.min())); del t,p,dx,vb,r
print(f'CENSUS-FINAL {sys.argv[1]}: below-gate {nb:,} worst {worst:.4f}')
assert nb==0, f'{nb} cells still below the 1 Hz p5 gate'
PYEOF

echo "=== 5/5 structure checks + promote ===" >> $LOG
$PY -u code/check_puml_faces.py $GATED >> $LOG 2>&1 || { echo ABORT-faces; exit 1; }
cp $GATED $FINAL
echo "ALT-FINALIZE-DONE $FINAL" | tee -a $LOG
