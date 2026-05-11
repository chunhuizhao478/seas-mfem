#!/bin/bash
#SBATCH -J bp5_1000m            # Job name
#SBATCH -o bp5_1000m_%j.out     # stdout (%j = job ID)
#SBATCH -e bp5_1000m_%j.err     # stderr
#SBATCH -p normal               # Queue
#SBATCH -N 4                    # nodes
#SBATCH -n 224                  # MPI tasks (56 per node)
#SBATCH -t 48:00:00             # wall time
#SBATCH -A EAR20006             # TACC allocation

export LC_ALL=C
export LANG=C

module load hypre/2.31.0
module load mumps/5.3
module load parmetis
# R-309 / plan §Phase 2b "Files to Modify": parallel HDF5 module is
# required by VTKHDF fault output (the production default on
# MFEM_USE_HDF5=YES builds).  Module name on Frontera follows the
# `phdf5/<version>` convention; some clusters prefer `hdf5/1.14.x`.
# Tolerate both via `||`.
module load phdf5 2>/dev/null || module load hdf5/1.14.4 2>/dev/null || true

# R-309 / plan §Phase 2d.3: opt-in ZFP lossy compression flags.
# Defaults (CLAUDE.md "ZFP lossy output"):
#   bulk:  1e-3 (m/s, velocity)
#   fault: 1e-12 (slip-rate floor)
# Uncomment AFTER H5Z-ZFP is verified on this cluster (Phase 2d.4 gate).
# export HDF5_PLUGIN_PATH=$WORK/h5z-zfp/lib/plugin
# PARAVIEW_BULK_ZFP="--paraview-bulk-zfp-tol 1e-3"
# PARAVIEW_FAULT_ZFP="--paraview-fault-zfp-tol 1e-12"
module list

# Ensure hypre/2.31.0 is found before PETSc's old bundled libHYPRE.so
export LD_LIBRARY_PATH=${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH}

cd /scratch2/10024/zhaochun/seas-project/seas-mfem/miniapps/seas

# Rebuild with latest optimizations
make seas_bp5_full

# Generate mesh if not present
if [ ! -f bp5/mesh/bp5_1000m.msh ]; then
    gmsh -3 bp5/mesh/bp5.geo -o bp5/mesh/bp5_1000m.msh -setnumber h 1.0
fi

mkdir -p bp5/results_1000m

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_1000m \
      --checkpoint-interval 5000 \
      --mumps
