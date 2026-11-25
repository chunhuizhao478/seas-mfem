#!/bin/bash
#SBATCH -J v_ucerf3d_pcg_4n
#SBATCH -o v_ucerf3d_pcg_4n_%j.out
#SBATCH -e v_ucerf3d_pcg_4n_%j.err
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 03:00:00
#SBATCH -A EAR20006

# UCERF 3-D CPU PCG timing on 4 Vista gg nodes (576 ranks): v_ucerf3d_pcg.sbatch
# with PART=run, on the binary the 1-node job 1170873 built at this code, and
# the meshes of Horizon job 106454 (gmsh 4.13.1 on macOS; sha256 below).  All
# settings are in this file, so nothing is overridden at submit time:
#   sbatch --dependency=afterok:1170873 jobs/vista_ucerf3d/v_ucerf3d_pcg_4node.sbatch
#   (from miniapps/seas of the Vista 3-D clone)
export PART=run
export PILOT_SHA256=9a4481a184d4211bdf0efc2e3fecd41ee1ccb28975c21f6021bec9a61a6d63cf
export TARGET_SHA256=47502639117c0d68346a05f0390ee1e6ecb204a5e47fb072d51ea616dffa7c11
exec bash "${SLURM_SUBMIT_DIR}/jobs/vista_ucerf3d/v_ucerf3d_pcg.sbatch"
