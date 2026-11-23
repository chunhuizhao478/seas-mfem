#!/bin/bash
#SBATCH -J h_d2_build
#SBATCH -o h_d2_build_%j.out
#SBATCH -e h_d2_build_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -A EAR26024

# GPU dynamic-rupture plan D2 on Horizon: jobs/vista/d2_gh_build.sbatch (the
# CUDA build of seas_nitsche_dyn_driver and seas_test_nitsche_sipg_mf, the
# unit checks on 1 GPU and on 4 ranks with one GPU each, and the Nsight
# trace on 1 GPU) on one GB200 node, with the settings of
# jobs/horizon/h_dyn_env.sh.  The first run also builds HDF5 into HDF5_ROOT
# and MFEM from scratch (device code for compute_90, machine code sm_100).
#
#   sbatch jobs/horizon/h_d2_build.sbatch   (from miniapps/seas of the Horizon
#                                            clone; or h_d2_submit.sh)

set -uo pipefail
source "${SLURM_SUBMIT_DIR}/jobs/horizon/h_dyn_env.sh" || exit 1
exec bash "${SLURM_SUBMIT_DIR}/jobs/vista/d2_gh_build.sbatch"
