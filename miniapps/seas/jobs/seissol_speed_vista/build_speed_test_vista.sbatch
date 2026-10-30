#!/bin/bash
#SBATCH -J sst_build
#SBATCH -o sst_build_%j.out
#SBATCH -e sst_build_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 03:00:00
#SBATCH -A EAR24034
#SBATCH --exclude=i617-131,i617-101
# SeisSol and MFEM speed test on Vista (plan rev1, Chunhui Zhao 2026-10-30), step B: every binary and mesh
# the runs need, built in one job at one pushed commit of this repository.
#   1. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   2. MFEM: build_vista.sh (MFEM + seas_nitsche_dyn_driver + the six gate tests; writes
#      miniapps/seas/BUILD_COMMIT_VISTA), and the two 200 m gold meshes checked against
#      jobs/vista_profiling/meshes.sha256 (copied from $WORK/seas-mfem-benchmarks if absent)
#   3. SeisSol elastic o2, o3, o4 for the Grace CPUs, in $SEISSOL_VISTA/build-o<N>-elastic-m0-sve128.
#      Toolchain = the Vista CPU kit's environment file $SEISSOL_VISTA/kit/seissol-env-vista.sh, pinned
#      by its sha256 below (gcc 13.2.0, openmpi 5.0.5, cmake 3.31.5, phdf5 1.14.6, parallel netcdf 4.9.2,
#      eigen 3.4.0, yaml-cpp 0.8.0; Lua, 64-bit ParMETIS, ASAGI, easi from the kit's deps/), source
#      $SEISSOL_VISTA/src/SeisSol at the kit's pin a8e7f04f, and the kit's CMake flags (cmake_args_for in
#      kit/build_seissol_vista.sh) with ORDER = 2, 3, 4 and EQUATIONS = elastic.  The kit's GPU builds,
#      its viscoelastic configuration and vista-gh/ are not touched.
#   4. PUMGen v1.1.3 (SeisSol's mesh converter; Gmsh msh2 reader, no SCOREC, no NetCDF) in
#      $SEISSOL_VISTA/pumgen, with the same toolchain
#   5. the two meshes converted (pumgen -s msh2; TPV205's tags 101/103/105 first renamed 1/3/5 by
#      retag_msh2.py) to $SEISSOL_VISTA/meshes/<name>.puml.h5, with the tet
#      count checked against MFEM's and the boundary tag histogram printed; MESHES.sha256 there lists
#      the .msh and .puml.h5 sha256s the run jobs check
#   6. check_parity.py for the three problems with the kit's Python (numpy, PyYAML, tomllib)
# Submit from miniapps/seas/jobs/seissol_speed_vista/ in $WORK/seas-mfem-profiling at the pushed commit.
# Cost: 1 node x <= 3 h x 0.33 SU = <= 1.0 SU (estimate: MFEM 25 min, SeisSol about 10 min per order,
# PUMGen 5 min).
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd $ROOT
for i in 1 2 3 4; do
   git fetch -q origin && break
   [ $i = 4 ] && { echo "ERROR: git fetch failed"; exit 2; }
   sleep $(( 5 + RANDOM % 20 ))
done
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
T0=$(date +%s)
fail=0

echo "=== 2. MFEM ==="
JOBS=64 bash build_vista.sh > build_vista_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "build_vista.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -5 build_vista_$SLURM_JOB_ID.log
[ $rc -eq 0 ] || exit $rc
cd miniapps/seas
for m in tpv205/mesh/tpv2053d_200m.msh tpv104/mesh/tpv104_200m.msh; do
   [ -f $m ] || cp "$WORK/seas-mfem-benchmarks/miniapps/seas/$m" $m
   grep -q "^$(sha256sum $m | cut -d' ' -f1)  $m\$" jobs/vista_profiling/meshes.sha256 \
      && echo "$m: sha256 matches the gold mesh" \
      || { echo "ERROR: $m sha256 $(sha256sum $m | cut -d' ' -f1) is not the gold mesh"; exit 1; }
done
SEAS=$PWD
SST=$SEAS/jobs/seissol_speed_vista

echo "=== 3. SeisSol elastic o2, o3, o4 ==="
export SEISSOL_VISTA=$WORK/seissol
KIT_ENV=$SEISSOL_VISTA/kit/seissol-env-vista.sh
KIT_ENV_SHA256=5df20ce33e9ec4ef9e6deec905ad113a34d5bf934c6021ba2ae7f6f38e5f7209
SEISSOL_COMMIT=a8e7f04f45514feb5e9cada8f5d5bcc5e6aa61f6
PUMGEN_TAG=v1.1.3
PUMGEN_COMMIT=6518bc6269a847b1d947f095eb04ec73af43ede0
[ "$(sha256sum $KIT_ENV | cut -d' ' -f1)" = $KIT_ENV_SHA256 ] \
   || { echo "ERROR: $KIT_ENV sha256 $(sha256sum $KIT_ENV | cut -d' ' -f1), expected $KIT_ENV_SHA256"; exit 2; }
(
set +u
SEISSOL_ENV_VERBOSE=1 source $KIT_ENV
[ "$SEISSOL_ENV_OK" = 1 ] || { echo "ERROR: the kit environment self-check failed"; exit 2; }
SRC=$SEISSOL_VISTA/src/SeisSol
[ "$(git -C $SRC rev-parse HEAD)" = $SEISSOL_COMMIT ] \
   || { echo "ERROR: $SRC is at $(git -C $SRC rev-parse HEAD), expected $SEISSOL_COMMIT"; exit 2; }
[ -z "$(git -C $SRC status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: $SRC has local changes"; git -C $SRC status --short; exit 2; }
echo "SeisSol source $SRC: $(git -C $SRC describe --tags --always) $(git -C $SRC rev-parse HEAD)"
MAN=$SEISSOL_VISTA/BUILD_MANIFEST-speedtest.txt
{ echo "# SeisSol builds for the speed test, job $SLURM_JOB_ID on $(hostname), $(date)"
  echo "# seas-mfem commit $(git -C $ROOT rev-parse HEAD); kit env sha256 $KIT_ENV_SHA256"
  echo "# modules: $(module -t list 2>&1 | tr '\n' ' ')"; } > $MAN
ok=1
for O in 2 3 4; do
   B=$SEISSOL_VISTA/build-o$O-elastic-m0-sve128
   mkdir -p $B
   s=$(date +%s)
   # the kit's cmake_args_for (kit/build_seissol_vista.sh) with ORDER=$O, EQUATIONS=elastic, no mechanisms
   ( cd $B && cmake -DCMAKE_BUILD_TYPE=Release -DORDER=$O -DPRECISION=double \
        -DEQUATIONS=elastic -DNUMBER_OF_MECHANISMS=0 \
        -DHOST_ARCH=sve128 -DGEMM_TOOLS_LIST=PSpaMM -DGRAPH_PARTITIONING_LIBS=parmetis \
        -DASAGI=ON -DNETCDF=ON -DNUMA_AWARE_PINNING=ON \
        -DPLASTICITY_METHOD=nb -DDR_QUAD_RULE=stroud -DLOG_LEVEL_MASTER=info \
        -DTESTING=OFF -DTESTING_PYTHON=OFF -DPython3_EXECUTABLE=$VENV/bin/python3 \
        -DCMAKE_C_FLAGS=-mcpu=neoverse-v2 -DCMAKE_CXX_FLAGS=-mcpu=neoverse-v2 \
        $(cat $SEISSOL_VISTA/.seissol-extra-cmake-args 2>/dev/null) $SRC ) > $B/configure_$SLURM_JOB_ID.out 2>&1 \
      && ( cd $B && make -j 128 ) > $B/make_$SLURM_JOB_ID.out 2>&1
   rc=$?
   BIN=$(find $B -maxdepth 1 -type f -perm -u+x -name 'seissol-elastic-o*' ! -name 'proxy*' | head -1)
   echo "o$O: exit $rc in $(( $(date +%s) - s )) s, binary ${BIN:-none}"
   if [ $rc -ne 0 ] || [ -z "$BIN" ]; then
      tail -30 $B/make_$SLURM_JOB_ID.out 2>/dev/null || tail -30 $B/configure_$SLURM_JOB_ID.out
      echo "o$O: FAILED (exit $rc)" >> $MAN; ok=0; continue
   fi
   [ "$(grep -m1 '^ORDER:' $B/CMakeCache.txt | cut -d= -f2)" = $O ] || { echo "ERROR: $B cache ORDER is not $O"; ok=0; }
   grep -q '^GEMM_TOOLS_LIST:.*=PSpaMM$' $B/CMakeCache.txt || { echo "ERROR: $B GEMM_TOOLS_LIST is not PSpaMM"; ok=0; }
   { echo "o$O: $BIN"; echo "   sha256 $(sha256sum $BIN | cut -d' ' -f1)"; } >> $MAN
done
echo "=== 4. PUMGen $PUMGEN_TAG ==="
P=$SEISSOL_VISTA/src/PUMGen
if [ ! -d $P/.git ]; then
   for i in 1 2 3 4; do
      git clone -q --branch $PUMGEN_TAG https://github.com/SeisSol/PUMGen.git $P && break
      rm -rf $P; [ $i = 4 ] && { echo "ERROR: git clone PUMGen failed"; exit 2; }
      sleep $(( 2 ** i ))
   done
   git -C $P submodule update --init --recursive || { echo "ERROR: PUMGen submodules"; exit 2; }
fi
[ "$(git -C $P rev-parse HEAD)" = $PUMGEN_COMMIT ] \
   || { echo "ERROR: $P is at $(git -C $P rev-parse HEAD), expected $PUMGEN_TAG $PUMGEN_COMMIT"; exit 2; }
PB=$P/build; rm -rf $PB; mkdir -p $PB
s=$(date +%s)
( cd $PB && cmake -DCMAKE_BUILD_TYPE=Release -DSCOREC=OFF -DSIMMETRIX=OFF -DNETCDF=OFF -DPARASOLID=OFF \
     -DCMAKE_INSTALL_PREFIX=$SEISSOL_VISTA/pumgen .. && make -j 64 && make install ) > $PB/build_$SLURM_JOB_ID.out 2>&1
rc=$?; echo "PUMGen build exit $rc in $(( $(date +%s) - s )) s"
[ $rc -eq 0 ] && [ -x $SEISSOL_VISTA/pumgen/bin/pumgen ] || { tail -30 $PB/build_$SLURM_JOB_ID.out; exit 2; }
echo "pumgen $PUMGEN_TAG: $SEISSOL_VISTA/pumgen/bin/pumgen sha256 $(sha256sum $SEISSOL_VISTA/pumgen/bin/pumgen | cut -d' ' -f1)" >> $MAN
echo "=== 5. meshes ==="
M=$SEISSOL_VISTA/meshes; mkdir -p $M
: > $M/MESHES.sha256.new
for spec in "tpv205/mesh/tpv2053d_200m.msh|tpv2053d_200m|1894583|101:1 103:3 105:5" "tpv104/mesh/tpv104_200m.msh|tpv104_200m|2466549|"; do
   IFS='|' read -r msh name ntet retag <<< "$spec"
   in=$SEAS/$msh
   if [ -n "$retag" ]; then
      # SeisSol reads a face tag 3 or >= 65 as a rupture face, so the TPV205 mesh's 101 (free surface) and
      # 105 (absorbing) are renamed 1 and 5 (and 103 -> 3) in a copy; see retag_msh2.py
      in=$M/${name}_seissol_tags.msh
      $VENV/bin/python3 -I $SST/retag_msh2.py $SEAS/$msh $in $retag || { ok=0; continue; }
   fi
   s=$(date +%s)
   ( cd $M && mpirun -np 1 $SEISSOL_VISTA/pumgen/bin/pumgen -s msh2 $in $name ) > $M/$name.pumgen_$SLURM_JOB_ID.log 2>&1
   rc=$?; echo "$name: pumgen exit $rc in $(( $(date +%s) - s )) s"
   [ $rc -eq 0 ] && [ -f $M/$name.puml.h5 ] || { tail -20 $M/$name.pumgen_$SLURM_JOB_ID.log; ok=0; continue; }
   n=$(h5ls $M/$name.puml.h5/connect | sed -n 's/.*Dataset {\([0-9]*\),.*/\1/p')
   [ "$n" = $ntet ] && echo "$name: $n tets = MFEM's" || { echo "ERROR: $name has '$n' tets, MFEM $ntet"; ok=0; }
   h5dump -d /boundary -b LE -o $M/$name.boundary.bin $M/$name.puml.h5 > /dev/null \
      && $VENV/bin/python3 -I -c "
import numpy as np, sys
b = np.fromfile(sys.argv[1], dtype='<i4')
t = np.concatenate([(b >> (8 * k)) & 0xff for k in range(4)])
u, c = np.unique(t[t > 0], return_counts=True)
print('   boundary face tags (tag: faces, fault faces counted from both sides):', ', '.join('%d: %d' % x for x in zip(u, c)))" $M/$name.boundary.bin
   rm -f $M/$name.boundary.bin ${M}/${name}_seissol_tags.msh
   { echo "$(sha256sum $SEAS/$msh | cut -d' ' -f1)  $msh"; echo "$(sha256sum $M/$name.puml.h5 | cut -d' ' -f1)  $name.puml.h5"; } >> $M/MESHES.sha256.new
done
mv $M/MESHES.sha256.new $M/MESHES.sha256
cat $M/MESHES.sha256 >> $MAN
echo "=== 6. deck parity (MFEM deck against SeisSol deck, all three problems) ==="
for P in tpv205 tpv104 tpv102; do
   $VENV/bin/python3 -I $SST/check_parity.py --problem $P --repo $ROOT > $M/parity_$P.txt 2>&1 \
      && tail -1 $M/parity_$P.txt || { cat $M/parity_$P.txt; echo "ERROR: $P decks differ (check_parity.py)"; ok=0; }
done
echo "=== manifest $MAN ==="; cat $MAN
[ $ok = 1 ]
) || fail=1
echo "total $(( $(date +%s) - T0 )) s"
[ $fail -eq 0 ] || { echo "ERROR: a step failed; do not submit the layout probe or the runs"; exit 1; }
echo "all builds and meshes done"
