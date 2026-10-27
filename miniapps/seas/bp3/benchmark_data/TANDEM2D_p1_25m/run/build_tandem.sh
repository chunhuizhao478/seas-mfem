#!/bin/bash
# Out-of-source Tandem builds (2D, p = 1 and 2) in new untracked folders of the existing
# checkout /Users/chunhuizhao/projects/tandem (a26753d6), Release, the tandem app only.
# Deps: tandem-petsc321 env (clone of mumps-petsc; mpicc/mpicxx = clang 22 + Open MPI 5.0.10,
# PETSc 3.21.6 real double with MUMPS + ParMETIS, zlib, OpenMP), Eigen3 3 + cmake 4.2.1 from mfem-dev,
# Lua 5.4.8 + pkg-config 0.29.2 from easi-harness, /usr/bin/python3 3.9.6 + numpy 2.0.2.
# ARCH=noarch: the Frontera builds' x86 ARCH values (hsw, skx, ...) do not exist on arm64.
# Holds the Mac run lock for the whole build.
set -o pipefail
source /Users/chunhuizhao/miniforge/etc/profile.d/conda.sh; conda activate tandem-petsc321
export OMP_NUM_THREADS=1
P=$CONDA_PREFIX
DEV=/Users/chunhuizhao/miniforge/envs/mfem-dev
EH=/Users/chunhuizhao/miniforge/envs/easi-harness
B=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/tandem-build
T=/Users/chunhuizhao/projects/tandem
LOCK=/Users/chunhuizhao/projects/.seas-mac-run.lock
have=0
release() { [ "$have" = 1 ] && grep -q "pid $$ " "$LOCK/owner" 2>/dev/null && rm -f "$LOCK/owner" && rmdir "$LOCK" && echo "lock released"; }
trap release EXIT
until mkdir "$LOCK" 2>/dev/null; do echo "$(date) waiting for lock: $(cat $LOCK/owner 2>/dev/null)"; sleep 60; done
have=1; echo "pid $$ thread 2d-inplane (Tandem 2D build) since $(date)" > "$LOCK/owner"
for p in ${DEGREES:-1 2}; do
  D=$T/build_2d_p$p; mkdir -p $D; cd $D
  t0=$(date +%s)
  echo "##### p=$p configure $(date)"
  set -x
  PKG_CONFIG_PATH=$P/lib/pkgconfig $DEV/bin/cmake $T \
    -DCMAKE_BUILD_TYPE=Release -DDOMAIN_DIMENSION=2 -DPOLYNOMIAL_DEGREE=$p -DARCH=noarch \
    -DCMAKE_C_COMPILER=$P/bin/mpicc -DCMAKE_CXX_COMPILER=$P/bin/mpicxx \
    -DPKG_CONFIG_EXECUTABLE=$EH/bin/pkg-config \
    -DEigen3_DIR=$DEV/share/eigen3/cmake \
    -DLUA_INCLUDE_DIR=$EH/include -DLUA_LIBRARY=$EH/lib/liblua.5.4.dylib \
    -DPython3_EXECUTABLE=/usr/bin/python3 \
    -DCMAKE_PREFIX_PATH=$P -DOpenMP_ROOT=$P > $B/configure_p$p.log 2>&1
  rc=$?; set +x
  echo "configure rc $rc ($(( $(date +%s) - t0 )) s)"; [ $rc -eq 0 ] || { tail -25 $B/configure_p$p.log; exit 1; }
  t1=$(date +%s)
  make -j6 tandem > $B/make_p$p.log 2>&1; rc=$?
  for x in "$D/app/tandem"; do echo "binary $x $(ls -l $x | awk "{print \$6, \$7, \$8}")"; otool -L $x | grep -E "not found|petsc|lua|parmetis|mpi\." ; $x --help > $B/help_p$p.txt 2>&1; echo "--help rc $? ($(wc -l < $B/help_p$p.txt) lines)"; done
  echo "make tandem rc $rc ($(( $(date +%s) - t1 )) s)"; [ $rc -eq 0 ] || { grep -E 'error|Error' $B/make_p$p.log | head -20; exit 1; }
done
