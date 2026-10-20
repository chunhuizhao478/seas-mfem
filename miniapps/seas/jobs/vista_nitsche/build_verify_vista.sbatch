#!/bin/bash
#SBATCH -J nitsche_build_verify
#SBATCH -o nitsche_build_verify_%j.out
#SBATCH -e nitsche_build_verify_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 03:00:00
#SBATCH -A EAR24034
# Vista build + verification for the Nitsche TPV205 / TPV104 benchmarks (Chunhui Zhao,
# 2026-10-16), with the TPV6 / TPV7 / TPV31 additions of plan rev 1 (2026-10-19): their
# meshes, the setup unit test, two 1000 m driver smokes and two more load checks.
# Network work on the login node, compiling here:
#   cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
#       https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
#   cd seas-mfem-benchmarks/miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch
# Steps:
#   0. the repo CLAUDE.md guard (GitHub first), plus no untracked sources
#   1. build_vista.sh: MFEM + seas_nitsche_dyn_driver + six gate tests (modules only)
#   2. meshes from the committed .geo with TACC gmsh (no OpenCASCADE, so every .geo uses the built-in
#      kernel), msh 2.2, sha256 recorded; the job stops if gmsh fails or a mesh has no elements:
#      tpv205/mesh/tpv2053d_1000m.msh, tpv2053d_200m.msh, tpv104/mesh/tpv104_200m.msh
#      (tpv104_200m.geo equals tpv102_200m.geo apart from comments), and
#      tpv6/mesh/tpv6_nitsche_1000m.msh, tpv6_nitsche_100m.msh (TPV6 and TPV7 share it),
#      tpv31/mesh/tpv31_nitsche_1000m.msh, tpv31_nitsche_200m.msh
#   3. gate tests, as on the Mac at c906a8a: gates 74/0, pf 22/0, hook-spy 16/16 at np 1, 2, 4 + abort,
#      lsw 13/13, srw 9/9, lsw-inline 14/14 at np 1, 2, 4; tpv-setup 30/30 at np 1 and 2 (cloud)
#   4. driver smoke: TPV205 1000 m p=1 to 1.2 s on 8 ranks.  Mac reference (same deck, its own
#      1000 m mesh): "816 steps, t = 1.20137 s; max|V| 1.52293 m/s; max|V|dt/Dc 0.0056054".  The
#      Vista mesh can differ from the Mac one (gmsh output depends on platform), so compare the
#      numbers, not bytes.  Then TPV6 1000 m p=1 to 6 s and TPV31 1000 m p=1 to 5 s, 8 ranks, with
#      the cloud numbers (4 ranks, gmsh 4.15.2 meshes) printed next to them.
#   5. 200 m load check: TPV205 then TPV104 200 m, p = 2, all 144 ranks, tfinal 0.003 s; node memory
#      peak, tets, implied min 6V/A and dt, steps to 12 s.  Confirms the ranks per node of the runs.
#      Then TPV31 200 m p = 2 on 144 ranks and TPV6 100 m p = 1 on 48 ranks, the ranks per node of the
#      TPV6/TPV7 100 m jobs (every rank reads the whole serial mesh, 3.2-3.4 GB per rank at this
#      size, so 144 per node would not fit).  On one node these two check only the driver's setup
#      (mesh read, partition, operator): after it the whole problem sits on this node, which the runs
#      spread over 4 to 48 nodes, and both ran out of memory there in job 1083417.  Memory on the
#      run layout: loadcheck_tpv6_100m_p1.sbatch and loadcheck_tpv31_200m_p1.sbatch.
# Exit status: 3 when a mesh fails (nothing after it would mean anything), otherwise 1 when a gate
# test or a smoke run fails or the TPV31/TPV6 load check dies before the driver's setup has ended
# (no dt line; all are still run and printed), else 0.  The TPV205 and TPV104 p = 2 load checks at
# 144 ranks run out of memory on one node (they did in every build so far) and do not count.
# Cost: 1 node x <= 3 h x 0.33 SU = <= 1 SU (the TPV6 100 m mesh took 306 s of gmsh in the cloud).
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd $ROOT
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
T0=$(date +%s)
JOBS=64 bash build_vista.sh > build_vista_$SLURM_JOB_ID.log 2>&1
rc=$?; echo "build_vista.sh exit $rc after $(( $(date +%s) - T0 )) s"; tail -20 build_vista_$SLURM_JOB_ID.log
[ $rc -eq 0 ] || exit $rc
# same module environment for the runs below
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
cd miniapps/seas
echo "=== meshes (gmsh $(gmsh --version 2>&1)) ==="
for g in tpv205/mesh/tpv2053d_1000m tpv205/mesh/tpv2053d_200m tpv104/mesh/tpv104_200m \
         tpv6/mesh/tpv6_nitsche_1000m tpv31/mesh/tpv31_nitsche_1000m tpv31/mesh/tpv31_nitsche_200m \
         tpv6/mesh/tpv6_nitsche_100m; do
   t=$(date +%s); gmsh -3 $g.geo -format msh22 -o $g.msh > $(basename $g)_gmsh_$SLURM_JOB_ID.log 2>&1
   rc=$?; ne=$(grep -a -A1 '^\$Elements' $g.msh 2>/dev/null | tail -1)
   echo "$g.msh: gmsh exit $rc in $(( $(date +%s) - t )) s, ${ne:-0} elements"
   [ $rc -eq 0 ] && [ "${ne:-0}" -gt 0 ] \
      || { echo "ERROR: $g.msh not made; $(basename $g)_gmsh_$SLURM_JOB_ID.log:"; grep -a "Error" $(basename $g)_gmsh_$SLURM_JOB_ID.log | head -5; exit 3; }
done
FAIL=0
sha256sum tpv205/mesh/tpv2053d_1000m.msh tpv205/mesh/tpv2053d_200m.msh tpv104/mesh/tpv104_200m.msh \
          tpv6/mesh/tpv6_nitsche_1000m.msh tpv31/mesh/tpv31_nitsche_1000m.msh tpv31/mesh/tpv31_nitsche_200m.msh \
          tpv6/mesh/tpv6_nitsche_100m.msh | tee MESH_SHA256_VISTA
echo "=== gate tests ==="
for t in test-nitsche-gates test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw test-nitsche-lsw-inline \
         test-nitsche-tpv-setup; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   nf=$(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $nf"
   [ $rc -eq 0 ] && [ $nf -eq 0 ] || FAIL=$((FAIL + 1))
done
echo "=== driver smoke: TPV205 1000 m p=1 1.2 s, 8 ranks ==="
OUT=$SCRATCH/nitsche_benchmarks/verify_tpv205_1000m_p1_1p2s_$SLURM_JOB_ID; mkdir -p $OUT
s=$(date +%s)
ibrun -n 8 -o 0 ./seas_nitsche_dyn_driver --config tpv205/configs/tpv205_nitsche.toml \
   --mesh tpv205/mesh/tpv2053d_1000m.msh --order 1 --tfinal 1.2 --output-dir $OUT --log-every 1000 \
   > $OUT/driver.log 2>&1
echo "driver exit $? in $(( $(date +%s) - s )) s"; grep "=== done" $OUT/driver.log || FAIL=$((FAIL + 1))
echo "Mac reference: 816 steps, t = 1.20137 s; max|V| 1.52293 m/s; max|V|dt/Dc 0.0056054"
for case in "tpv6 tpv6/configs/tpv6_nitsche.toml tpv6/mesh/tpv6_nitsche_1000m.msh 6.0" \
            "tpv31 tpv31/configs/tpv31_nitsche.toml tpv31/mesh/tpv31_nitsche_1000m.msh 5.0"; do
   set -- $case
   echo "=== driver smoke: $1 1000 m p=1 $4 s, 8 ranks ==="
   OUT=$SCRATCH/nitsche_benchmarks/verify_$1_1000m_p1_$SLURM_JOB_ID; mkdir -p $OUT
   s=$(date +%s)
   ibrun -n 8 -o 0 ./seas_nitsche_dyn_driver --config $2 --mesh $3 --order 1 --tfinal $4 \
      --output-dir $OUT --log-every 1000 > $OUT/driver.log 2>&1
   echo "driver exit $? in $(( $(date +%s) - s )) s"; grep "=== done" $OUT/driver.log || FAIL=$((FAIL + 1))
done
echo "cloud reference TPV6 (4 ranks, built-in kernel mesh): 4081 steps, t = 6.00044 s; max|V| 9.29284 m/s; max|V|dt/Dc 0.034159"
echo "cloud reference TPV31 (4 ranks, built-in kernel mesh): 3375 steps, t = 5.00027 s; max|V| 1.5991 m/s; max|V|dt/Dc 0.0131621"
# Load checks: every rank reads the whole serial mesh at startup
# (drivers/nitsche_dyn_driver.cpp:829), so 144 ranks on one node may not fit in its memory.
# tfinal 0.003 s (about 25 steps at p = 2).  Node memory in use (MemTotal - MemAvailable,
# /proc/meminfo) is sampled every 2 s.  The driver prints no h_min, so min h_e / c_p,e is
# derived from its dt line: dt p^2 / (cfl C(p)), cfl 0.5; times 6000 m/s it is h_min for
# TPV205 and TPV104 (c_p = 6000 m/s everywhere; 26.81 m for the TPV102 200 m mesh).  TPV6 and
# TPV31 have more than one c_p, so only min h/c_p counts there.
# TPV205 (smaller mesh) first, so its numbers survive if TPV104 runs out of memory.
for case in "tpv205 tpv205/configs/tpv205_nitsche.toml tpv205/mesh/tpv2053d_200m.msh 2 144" \
            "tpv104 tpv104/configs/tpv104_nitsche.toml tpv104/mesh/tpv104_200m.msh 2 144" \
            "tpv31 tpv31/configs/tpv31_nitsche.toml tpv31/mesh/tpv31_nitsche_200m.msh 2 144" \
            "tpv6 tpv6/configs/tpv6_nitsche.toml tpv6/mesh/tpv6_nitsche_100m.msh 1 48"; do
   set -- $case; PROB=$1; DECK=$2; MESH=$3; P=$4; NR=$5
   echo "=== load check: $PROB $MESH p=$P, $NR ranks, tfinal 0.003 s ==="
   OUT=$SCRATCH/nitsche_benchmarks/loadcheck_${PROB}_p${P}_$SLURM_JOB_ID; mkdir -p $OUT
   ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
   MON=$!
   s=$(date +%s)
   ibrun -n $NR -o 0 ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH --order $P --tfinal 0.003 \
      --output-dir $OUT --log-every 1 > $OUT/driver.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   echo "driver exit $rc in $(( $(date +%s) - s )) s"
   # the driver prints its dt line once the setup has ended (drivers/nitsche_dyn_driver.cpp:1668)
   [ $rc -eq 0 ] || [ $PROB = tpv205 ] || [ $PROB = tpv104 ] || grep -q -a "  dt = " $OUT/driver.log \
      || FAIL=$((FAIL + 1))
   grep -a "tets;" $OUT/driver.log | head -1
   grep -a "  dt = " $OUT/driver.log | head -1
   grep -a "  dt = " $OUT/driver.log | head -1 | awk '{for(i=1;i<=NF;i++){if($i=="dt"&&$(i+1)=="=")dt=$(i+2); if($i=="C(p)")C=$(i+1)}} END{if(dt==""||C==""){print "no dt line in driver.log"; exit} r=dt*p*p/(0.5*C); printf "implied min h/c_p %.4e s (x 6000 m/s = %.2f m); steps to 12 s %d\n", r, r*6000, int(12/dt+0.999)}' p=$P
   grep -a "=== done" $OUT/driver.log || tail -15 $OUT/driver.log
   awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "node memory: peak in use %.1f GB of %.1f GB (MemTotal - MemAvailable, 2 s samples)\n", m, t}' $OUT/meminfo_kB.log
done
sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,Elapsed,MaxRSS,MaxVMSize,NTasks 2>&1 | head -20
echo "(MaxRSS is per task and appears only for srun steps; ibrun with Open MPI may leave it empty, so the /proc/meminfo peak above is the number that decides the ranks per node)"
echo "total $(( $(date +%s) - T0 )) s"
echo "failed checks: $FAIL"
[ $FAIL -eq 0 ] || exit 1
