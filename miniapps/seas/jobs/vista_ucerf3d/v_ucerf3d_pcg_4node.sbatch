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
# with PART=all: it rebuilds seas_test_qd_pcg at this commit (PETSc and MFEM
# reused when unchanged), so the binary and the timings carry the job's own
# commit, then runs on the meshes of Horizon job 106454 (gmsh 4.13.1 on macOS;
# sha256 below).  PART=run cannot follow a pull: the commit stamp
# cycle/.commit.cache is a make prerequisite, so make -q calls the binary
# stale (job 1170879, exit 2).  All settings are in this file, so nothing is
# overridden at submit time:
#   sbatch jobs/vista_ucerf3d/v_ucerf3d_pcg_4node.sbatch
#   (from miniapps/seas of the Vista 3-D clone)
export PART=all
export PILOT_SHA256=9a4481a184d4211bdf0efc2e3fecd41ee1ccb28975c21f6021bec9a61a6d63cf
export TARGET_SHA256=47502639117c0d68346a05f0390ee1e6ecb204a5e47fb072d51ea616dffa7c11
exec bash "${SLURM_SUBMIT_DIR}/jobs/vista_ucerf3d/v_ucerf3d_pcg.sbatch"
