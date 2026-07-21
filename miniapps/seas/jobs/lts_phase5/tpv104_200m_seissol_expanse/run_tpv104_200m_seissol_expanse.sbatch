#!/bin/bash
# =============================================================================
# LTS PHASE 5 — SeisSol side of the TPV104-200m head-to-head  --  SDSC Expanse
#
# SeisSol v1.1.3 o4 (= p3 polynomials, order-matched to the MFEM p3/ADER-O4 run)
# on the SAME 200 m mesh as the MFEM Phase-5 run: tpv104.puml.h5 is pumgen's
# conversion of a byte-for-byte copy of miniapps/seas/tpv104/mesh/tpv104_200m.msh
# (2,464,689 tets, 37,904 DR faces).  EndTime = 2.0 s, ClusteredLTS = 2, pure
# Godunov flux — the partner numbers are the MFEM GTS (~0.27 sim-s/h) and MFEM
# LTS (jobs 52332649 / 52340525) runs from ../tpv104_200m_lts_speed_expanse/.
#
# SINGLE-SHOT: no checkpoint chain (2 s window; Checkpoint off in the .par).
# Setup follows run_safs_v4_0_0_pref_case1_small_o4.sbatch (the SAFS SeisSol
# deck family): same account/binary/env block (every var is a hard-won fix; do
# not trim), same outputs/ + STOP_RUN/.SIMULATION_DONE conventions.  The chain/
# checkpoint machinery is intentionally absent (2 s window, Checkpoint=0).
#
# DEPLOY (rsync from the laptop, like the other SeisSol decks — NOT run from
# the seas-mfem repo tree on Expanse):
#   rsync -avP <repo>/miniapps/seas/jobs/lts_phase5/tpv104_200m_seissol_expanse/ \
#         czhao1@login.expanse.sdsc.edu:/expanse/lustre/projects/lbl107/czhao1/safs/tpv104_200m_seissol/
#   rsync -avP ~/Downloads/tpv104/tpv104.puml.h5 \
#         czhao1@login.expanse.sdsc.edu:/expanse/lustre/projects/lbl107/czhao1/safs/tpv104_200m_seissol/
#   cd /expanse/lustre/projects/lbl107/czhao1/safs/tpv104_200m_seissol
#   sbatch run_tpv104_200m_seissol_expanse.sbatch      # or: sbatch -A lbl107 ...
#
# WATCH:  squeue -u $USER ; tail -f tpv104-ss-o4.<jobid>.out
# =============================================================================
#SBATCH --job-name=tpv104-ss-o4
#SBATCH --account=usc143             # the account the SeisSol qstore build ran under; override: sbatch -A lbl107
#SBATCH --partition=compute
#SBATCH --constraint="lustre"        # REQUIRED: this run lives on /expanse/lustre
#SBATCH --nodes=2                    # SAME node count as the MFEM Phase-5 run (fair per-node throughput)
#SBATCH --ntasks-per-node=8          # 8 ranks x 16 cores = 128 cores/node = THE WHOLE NODE.
                                     # EPYC 7742 = 8 NUMA x 16 cores -> 1 rank per NUMA domain.
#SBATCH --cpus-per-task=16
#SBATCH --mem=140000M
#SBATCH -t 06:00:00                  # 2.0 sim-s on 2.46M tets at o4+LTS; generous margin
#SBATCH -o tpv104-ss-o4.%j.out
#SBATCH -e tpv104-ss-o4.%j.out
#SBATCH --export=ALL

set -o pipefail

# never dump core: SeisSol's LTS-watchdog abort() writes ~2.5 GB/rank otherwise.
ulimit -c 0

SEISSOL_ROOT=/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build
# o4 = ORDER 4 = p3.  Matches the MFEM side's order=3 + ader_order=4 exactly.
# The speed comparison is VOID on any other order — do not swap in o5/o6.
SEISSOL_EXE=${SEISSOL_EXE:-$SEISSOL_ROOT/SeisSol/build-rome-o4-elastic/seissol-elastic-o4-f64}

RUN_DIR=${SLURM_SUBMIT_DIR:-$PWD}
cd "$RUN_DIR" || exit 1
mkdir -p outputs                     # SeisSol will NOT create it; OutputFile='outputs/tpv104'

[ -x "$SEISSOL_EXE" ] || { echo "FATAL: no executable at $SEISSOL_EXE"; exit 2; }

# --- stop conditions (deck-family convention) --------------------------------
if [ -f outputs/.SIMULATION_DONE ]; then
  echo ">>> Simulation already reached EndTime. Nothing to do."; exit 0
fi
if [ -f STOP_RUN ]; then
  echo ">>> STOP_RUN present -- halted by user."; exit 0
fi

# --- deck pre-flights --------------------------------------------------------
# 1) the mesh must be staged (103 MiB — over GitHub's limit, so it travels by scp):
if [ ! -f tpv104.puml.h5 ]; then
    echo "FATAL: tpv104.puml.h5 missing.  scp it from the laptop:"
    echo "  scp /Users/chunhuizhao/Downloads/tpv104/tpv104.puml.h5 <expanse>:$RUN_DIR/"
    exit 2
fi
# 2) integrity pin (the conversion of the exact MFEM 200m mesh):
MD5="$(md5sum tpv104.puml.h5 2>/dev/null | cut -d' ' -f1)"
if [ -n "$MD5" ] && [ "$MD5" != "5911b0b4ec29a8f61968d5b32ed0a11c" ]; then
    echo "FATAL: tpv104.puml.h5 md5 $MD5 != pinned 5911b0b4ec29a8f61968d5b32ed0a11c"
    echo "       (must be the conversion of the MFEM tpv104_200m.msh — see MESH provenance in README)"
    exit 2
fi
[ -z "$MD5" ] && echo "WARN: md5sum unavailable — mesh integrity pin SKIPPED"
# 3) parity guards — a drifted .par must not silently change the physics window.
#    ([^0-9]|$) terminators: reject drift like t_0=1.05 / EndTime=2.0e1.
grep -qE '^t_0 *= *1\.0([^0-9]|$)'     parameters.par || { echo "FATAL: parameters.par t_0 != 1.0 (parity with MFEM T_nuc)"; exit 2; }
grep -qE '^EndTime *= *2\.0([^0-9]|$)' parameters.par || { echo "FATAL: parameters.par EndTime != 2.0 (speed-run window)"; exit 2; }
#    strip '!' comments first — the deck's own comments MENTION the key (an
#    unstripped grep would false-positive and kill every submission):
sed 's/!.*//' parameters.par | grep -qi 'numfluxnearfault' && { echo "FATAL: numfluxnearfault key present — needs SeisSol >= v1.2.0, binary is v1.1.3"; exit 2; }
for f in tpv104_fault.yaml tpv104_material.yaml tpv104_receivers.dat tpv104_faultreceivers.dat; do
    [ -f "$f" ] || { echo "FATAL: missing $f"; exit 2; }
done
echo "=== deck pre-flights OK (mesh md5 pinned, t_0=1.0, EndTime=2.0, v1.1.3-safe par)"

echo "############ TPV104 200m SeisSol o4 speed run | job $SLURM_JOB_ID"
echo "=== nodes=$SLURM_NNODES ranks=$SLURM_NTASKS cpus/rank=$SLURM_CPUS_PER_TASK"
echo "=== exe: $SEISSOL_EXE"
echo "=== dir: $RUN_DIR"

# --- toolchain (MUST match the build) ----------------------------------------
module purge
module load slurm cpu/0.17.3b gcc/10.2.0
module load openmpi/4.1.3 || module load openmpi/4.1.1
module list 2>&1

# --- runtime env (verbatim from the proven SAFS decks) -----------------------
# HDF5 >=1.10 POSIX file locking FAILS on Lustre (errno 524) and kills PUML/XDMF IO.
export HDF5_USE_FILE_LOCKING=FALSE
# 1 core per rank is reserved for SeisSol's communication thread -> OMP gets N-1.
export OMP_NUM_THREADS=$((SLURM_CPUS_PER_TASK - 1))
export OMP_PLACES=cores
export OMP_PROC_BIND=close
# ASYNC_MODE=SYNC (not THREAD): vtkhdf surface writer on the async I/O thread collided
# with the main-thread HDF5 receiver writer inside non-threadsafe parallel HDF5 1.14
# ('H5VL_reset_vol_wrapper' abort) at the first surface/receiver sync coincidence.
export ASYNC_MODE=SYNC
export ASYNC_BUFFER_ALIGNMENT=8388608
export XDMFWRITER_ALIGNMENT=8388608
export XDMFWRITER_BLOCK_SIZE=8388608
# ASAGI mode: no ASAGI grids in TPV104 (analytic material) — kept verbatim for
# env parity with the SAFS deck family (SEISSOL_ prefix required; harmless here).
export SEISSOL_ASAGI_MPI_MODE=OFF
# comm thread off: confirmed fix for the 64-rank LTS ghost timeout; harmless at 16 ranks.
export SEISSOL_COMMTHREAD=0

# --- UCX robustness (ride out fabric hiccups under MPI_THREAD_MULTIPLE) ------
export UCX_RC_MLX5_TIMEOUT=35000000.00us
export UCX_RC_MLX5_RNR_TIMEOUT=35000000.00us
export UCX_RC_MLX5_RETRY_COUNT=180
export UCX_RC_MLX5_RNR_RETRY_COUNT=180
export UCX_DC_MLX5_TIMEOUT=35000000.00us
export UCX_DC_MLX5_RNR_TIMEOUT=35000000.00us
export UCX_DC_MLX5_RETRY_COUNT=180
export UCX_DC_MLX5_RNR_RETRY_COUNT=180
export UCX_UD_MLX5_TIMEOUT=35000000.00us
export UCX_UD_MLX5_RETRY_COUNT=180

# *** SLURM 23.02 DOES NOT propagate --cpus-per-task from sbatch to srun. ***
# Without BOTH the env var and the explicit -c below, each rank gets ONE cpu and the
# run crawls at ~1/15 speed ('OpenMP worker affinity (this process): ( 0 )').
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

echo "=== OMP_NUM_THREADS=$OMP_NUM_THREADS  SRUN_CPUS_PER_TASK=$SRUN_CPUS_PER_TASK"
numactl -H 2>/dev/null | head -3

# --- run ---------------------------------------------------------------------
echo "=========================== SEISSOL START $(date) ==========================="
START=$SECONDS
# --unbuffered + stdbuf: SeisSol's logger never flushes under SLURM otherwise.
srun --mpi=pmix -n "$SLURM_NTASKS" -c "$SLURM_CPUS_PER_TASK" --cpu-bind=cores --unbuffered \
     stdbuf -oL -eL "$SEISSOL_EXE" parameters.par
rc=$?
ELAPSED=$((SECONDS - START))
echo "=========================== SEISSOL END   $(date)  exit=$rc  (${ELAPSED}s) ==========="

# --- speed metric + sanity readout -------------------------------------------
if [ $rc -eq 0 ]; then
    touch outputs/.SIMULATION_DONE     # deck-family convention: mark EndTime reached
    # sim-seconds per wall-hour over the whole run (init included), the MFEM-comparable number:
    SIMS_PER_HR=$(awk -v e="$ELAPSED" 'BEGIN { if (e > 0) printf "%.3f", 2.0 * 3600.0 / e }')
    echo "SPEED METRIC: 2.0 sim-s in ${ELAPSED}s wall -> ${SIMS_PER_HR} sim-s/wall-hour (init included)"
    echo "  Partners: MFEM GTS ~0.27 sim-s/h | MFEM LTS run 52340525 (see its .out)"
    echo "  For the init-free rate, subtract the wall time of the last 'init' log line."
else
    echo "SPEED METRIC: UNAVAILABLE (non-zero exit $rc) — check the log tail above."
fi
# v1.1.3 SILENTLY drops unlocated receivers — verify all 6 + 9 were found:
echo "=== receiver location check (expect 6 off-fault + 9 on-fault):"
ls -1 outputs/tpv104-receiver-* 2>/dev/null | wc -l
ls -1 outputs/tpv104-faultreceiver-* 2>/dev/null | wc -l
# LTS clustering the binary actually chose (for the report):
grep -iE "cluster|time step" "tpv104-ss-o4.${SLURM_JOB_ID}.out" 2>/dev/null | head -10 || true
echo "=== outputs:"; ls -lh outputs/ 2>/dev/null | head -15
exit $rc
