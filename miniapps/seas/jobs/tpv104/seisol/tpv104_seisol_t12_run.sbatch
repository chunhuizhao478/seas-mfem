#!/bin/bash
#SBATCH -J tpv104_seisol_t12
#SBATCH -o tpv104_seisol_t12_%j.out
#SBATCH -e tpv104_seisol_t12_%j.err
#SBATCH -p normal
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 12:00:00
#SBATCH -A EAR20006

# TPV104 SeisSol production run — t = 12 s, 500 m fault resolution,
# 400 ranks, full-domain mesh.
#
# Purpose
# -------
# Independent SeisSol reference run on the SAME geometric setup as
# MFEM's t12 dev run (results_t12_dev_job7677022).  Identical:
#   - Domain extent (60×60×60 km half-space)
#   - Fault geometry (vertical, y=0, 36×18 km sliding patch)
#   - Hypocenter location (0, 0, -7500 m)
#   - SCEC TPV104 parameter set (rs_f0, rs_b, rs_sr0, rs_iniSlipRate*,
#     RS_muW=0.2, t_0=0.5)  — these are the upstream SeisSol/Examples
#     reference parameters that generated benchmark_data/seisol/.
# Different:
#   - Mesh: full-domain unstructured tet (paired with MFEM's tpv104_500m)
#     vs upstream's half-mesh-and-mirror structured.
#   - SeisSol numerical scheme: ADER-DG, p=4, GPwise=1.
#
# Submission flow
# ---------------
# Pre-requisite: tpv104_seisol_mesh_build.sbatch must have produced
#   tpv104_seisol_500m.h5 + tpv104_seisol_500m.xdmf in the same
#   directory this script runs in.  Verified at startup.
#
# Cost
# ----
# Same wall-clock as MFEM's TPV104 production at 400 ranks: estimated
# ~6-10 h on Frontera normal-queue.  Allocated 12 h to be safe.
#
# After completion, the output station files at
#   output/tpv104_seisol-faultreceiver-*  (9 files)
#   output/tpv104_seisol-receiver-*       (6 files)
# can be diff'd directly against MFEM's tpv104_dthalf_station_*.dat
# via the visualizer at miniapps/seas/tpv104/visualize_results.py
# (will need a SeisSol-output adapter — see README.md).
#
# Toolchain (align with reference_seissol_frontera.md memory):
#   ICX 2023.1 + IMPI 2021.9 + GCC 9.1 + cmake 3.24
#   The `module load` block below uses Frontera's TACC modules; cross-
#   check against the working SeisSol Frontera build before submission.

export LC_ALL=C
export LANG=C

module purge
module load intel/23.1.0       2>/dev/null || module load intel/19.1.1
module load impi/21.9          2>/dev/null || module load impi/19.0.9
module load gcc/9.1
module load python3
module load hdf5/1.14.4        2>/dev/null || module load phdf5
module load cmake/3.24.2       2>/dev/null || module load cmake

# SeisSol install — adjust to match user's Frontera build location.
# The reference_seissol_frontera.md memory says SeisSol is built at
# /scratch2/.../seas-project/SeisSol — adjust if different.
SEISSOL_BIN="${SEISSOL_BIN:-/scratch2/10024/zhaochun/seas-project/SeisSol/build/SeisSol_Release_dskx_5_elastic}"
if [ ! -x "${SEISSOL_BIN}" ]; then
   echo "ERROR: SEISSOL_BIN=${SEISSOL_BIN} not executable."
   echo "       Set SEISSOL_BIN env var or rebuild SeisSol per"
   echo "       reference_seissol_frontera.md (ICX 2023.1 + IMPI 2021.9 + ORDER=5)."
   exit 1
fi

# Verify mesh build artefacts present.
for f in tpv104_seisol_500m.h5 tpv104_seisol_500m.xdmf \
         tpv104_seisol_material.yaml tpv104_seisol_fault.yaml \
         tpv104_seisol_faultreceivers.dat tpv104_seisol_receivers.dat \
         tpv104_seisol.par; do
   if [ ! -s "$f" ]; then
      echo "ERROR: required input file '$f' missing or empty in $(pwd)."
      echo "       Run tpv104_seisol_mesh_build.sbatch first; copy yaml +"
      echo "       receivers + par from miniapps/seas/jobs/tpv104/seisol/."
      exit 1
   fi
done

mkdir -p output

echo "========================================"
echo "SeisSol binary: ${SEISSOL_BIN}"
echo "Mesh: tpv104_seisol_500m  (PUMGen full-domain)"
echo "Target: tfinal = 12 s on 400 ranks"
echo "Reference parameters: RS_muW=0.2, t_0=0.5 (upstream Examples)"
echo "========================================"

ibrun "${SEISSOL_BIN}" tpv104_seisol.par

# Post-run sanity gate
RESULT_FILE="RESULT.txt"
NPP=$(ls -1 output/tpv104_seisol-faultreceiver-* 2>/dev/null | wc -l)
NREC=$(ls -1 output/tpv104_seisol-receiver-*    2>/dev/null | wc -l)

NAN_COUNT=0
for f in output/tpv104_seisol-faultreceiver-* output/tpv104_seisol-receiver-*; do
   if [ -s "$f" ] && grep -qiE "nan|inf" "$f"; then
      NAN_COUNT=$((NAN_COUNT + 1))
   fi
done

{
   echo "TPV104 SeisSol production run"
   echo "Mesh:                  tpv104_seisol_500m (full-domain)"
   echo "fault receivers found: ${NPP}"
   echo "off-fault receivers:   ${NREC}"
   echo "files with NaN/Inf:    ${NAN_COUNT}"
   if [ "${NPP}" -ge 9 ] && [ "${NAN_COUNT}" -eq 0 ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
