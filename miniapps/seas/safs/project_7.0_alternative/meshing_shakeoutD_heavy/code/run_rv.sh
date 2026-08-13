#!/bin/zsh
# run_rv.sh -- rv repair of the stage-1 mesh, detached.
#
# GUARD CHOICE: --cvm/--muscal (the 1 Hz gate, Vs/0.8), NOT --target-field.
# --target-field is the right bound only when the mesh was refined against that
# field. This mesh was adapted by mmg to the 1400-floored metric and is being
# shipped against the frequency gate, and h_target = min(h_fault, Vs/0.8, hmax)
# is <= Vs/0.8 EVERYWHERE -- so guarding on it would over-constrain apex moves
# and block repairs that are perfectly safe for what we actually deliver.
#
# --qual-tol 0.05 is per-tet, against each tet's OWN eta and inradius. Guarding
# against the mesh's GLOBAL worst eta instead let dt fall 10.8x (16.129 -> 1.487
# us) in an earlier campaign.
#
# Slivers here are 100 % fault-adjacent (97.4 % within 250 m, zero beyond 2 km),
# so this stage works entirely on tets with a frozen face: the three fault
# vertices cannot move and only the free apex can. Expect modest gains.
set -u
cd "${0:a:h}/.."
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
B=build_tmp
LOG=$B/rv.log
MESH=$B/shakeoutD_alt_heavy_s1.puml.h5
OUT=$B/shakeoutD_alt_heavy_s1_rv.puml.h5
CVM=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/seisol_quakeworx/safs_seisol_v3_2_0_RSSRW_ALTv2_freesurf/safs_material_cvm.nc
MUSCAL=/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc

[ -f $MESH ] || { echo "MISSING $MESH"; exit 1; }
if pgrep -f "leb_gate_close.py|relax_rv.py" >/dev/null 2>&1 || pgrep -x mmg3d_O3 >/dev/null 2>&1; then
  echo "REFUSING: a heavy job is already running."; exit 1
fi

echo "[start] $(date)"
nohup $PY -c "
import subprocess, sys
subprocess.Popen(
    [sys.argv[1], '-u', 'code/relax_rv.py',
     '--mesh', '$MESH', '--out', '$OUT',
     '--cvm', '$CVM', '--muscal', '$MUSCAL',
     '--gate', '0.8', '--qual-tol', '0.05'],
    stdout=open('$LOG', 'ab'), stderr=subprocess.STDOUT,
    start_new_session=True)
" $PY >> $LOG 2>&1 &!
sleep 4
pid=$(pgrep -f "code/relax_rv.py" | head -1)
[ -z "$pid" ] && { echo "[FAILED] nothing launched"; tail -15 $LOG; exit 1; }
echo "[launched] pid $pid -> $LOG"
