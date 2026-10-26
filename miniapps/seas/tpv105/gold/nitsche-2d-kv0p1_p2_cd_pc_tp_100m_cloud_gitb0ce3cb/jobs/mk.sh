#!/bin/bash
# cloud build wrapper for miniapps/seas (2d-inplane thread)
S=/tmp/claude-0/-home-user-seas-mfem/ce558891-cad7-5003-a322-13bbf23d3889/scratchpad
cd /home/user/seas-mfem/miniapps/seas
exec make -f cycle/cycle.mk MFEM_BUILD_DIR=$S/mfem-build SEAS_FS_LIB= SEAS_INCLUDES=-I/home/user/seas-mfem \
  HDF5_INCFLAGS="-I/usr/include/hdf5/openmpi -DH5Treclaim=H5Dvlen_reclaim" \
  HDF5_LIBS="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5" \
  "MFEM_MPIEXEC=mpirun --allow-run-as-root --oversubscribe" "$@"
