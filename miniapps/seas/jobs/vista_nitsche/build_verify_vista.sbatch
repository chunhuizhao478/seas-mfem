#!/bin/bash
#SBATCH -J nitsche_build_verify
#SBATCH -o nitsche_build_verify_%j.out
#SBATCH -e nitsche_build_verify_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
# Vista build + verification for the Nitsche TPV205 / TPV104 benchmarks (Chunhui Zhao,
# 2026-10-16).  Submit from the checkout's miniapps/seas/jobs/vista_nitsche/ after
#   cd $WORK && git clone https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
#   cd seas-mfem-benchmarks && git checkout <pushed nitsche-benchmarks commit>
# (network work on the login node; compiling here).  Steps:
#   1. build_vista.sh: MFEM + seas_nitsche_dyn_driver + six gate tests (modules only)
#   2. meshes from the committed .geo with TACC gmsh, msh 2.2, sha256 recorded:
#      tpv205/mesh/tpv2053d_1000m.msh, tpv2053d_200m.msh, tpv104/mesh/tpv104_200m.msh
#   3. gate tests (as on the Mac at c906a8a): gates 74/0, p5 14/0, pf 22/0, hook-spy 16/16 x3 + abort,
#      lsw 13/13, srw 9/9
#   4. driver smoke: TPV205 1000 m p=1 to 1.2 s on 8 ranks.  Mac reference (same deck, its own
#      1000 m mesh): "816 steps, t = 1.20137 s; max|V| 1.52293 m/s; max|V|dt/Dc 0.0056054".  The
#      Vista mesh can differ from the Mac one (gmsh output depends on platform), so compare the
#      numbers, not bytes.
#   5. 200 m load check: TPV104 and TPV205 200 m, p = 2, all 144 ranks, tfinal 0.003 s; node memory peak,
#      tets, implied min 6V/A and dt, steps to 12 s.  Decides 144 vs 128 ranks per node for the runs.
# Cost: 1 node x <= 2 h x 0.33 SU = <= 0.66 SU.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR/../../../.." || exit 2          # repository root
ROOT=$(pwd)
echo "repo $ROOT  commit $(git rev-parse HEAD)"
DIRTY=$(git status --porcelain --untracked-files=no | wc -l)
[ "$DIRTY" -eq 0 ] || { echo "ERROR: $DIRTY tracked files modified; run only from a clean pushed commit"; git status --short; exit 2; }
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
for g in tpv205/mesh/tpv2053d_1000m tpv205/mesh/tpv2053d_200m tpv104/mesh/tpv104_200m; do
   t=$(date +%s); gmsh -3 $g.geo -format msh22 -o $g.msh > $(basename $g)_gmsh_$SLURM_JOB_ID.log 2>&1
   echo "$g.msh: gmsh exit $? in $(( $(date +%s) - t )) s"
done
sha256sum tpv205/mesh/tpv2053d_1000m.msh tpv205/mesh/tpv2053d_200m.msh tpv104/mesh/tpv104_200m.msh | tee MESH_SHA256_VISTA
echo "=== gate tests ==="
for t in test-nitsche-gates test-nitsche-p5 test-nitsche-pf test-nitsche-hook-spy test-nitsche-lsw test-nitsche-srw; do
   s=$(date +%s); make MFEM_DIR=../.. HDF5_DIR="$TACC_PHDF5_DIR" SEAS_FS_LIB= $t > ${t}_$SLURM_JOB_ID.log 2>&1; rc=$?
   echo "$t exit $rc in $(( $(date +%s) - s )) s: PASSED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c PASSED), FAILED $(grep -v echo ${t}_$SLURM_JOB_ID.log | grep -c FAILED)"
done
echo "=== driver smoke: TPV205 1000 m p=1 1.2 s, 8 ranks ==="
OUT=$SCRATCH/nitsche_benchmarks/verify_tpv205_1000m_p1_1p2s_$SLURM_JOB_ID; mkdir -p $OUT
s=$(date +%s)
ibrun -n 8 -o 0 ./seas_nitsche_dyn_driver --config tpv205/configs/tpv205_nitsche.toml \
   --mesh tpv205/mesh/tpv2053d_1000m.msh --order 1 --tfinal 1.2 --output-dir $OUT --log-every 1000 \
   > $OUT/driver.log 2>&1
echo "driver exit $? in $(( $(date +%s) - s )) s"; grep "=== done" $OUT/driver.log
echo "Mac reference: 816 steps, t = 1.20137 s; max|V| 1.52293 m/s; max|V|dt/Dc 0.0056054"
# 200 m load check: every rank reads the whole serial mesh at startup
# (drivers/nitsche_dyn_driver.cpp:750), so 144 ranks on one node may not fit in its memory.
# p = 2, all 144 ranks, tfinal 0.003 s (about 25 steps).  Node memory in use
# (MemTotal - MemAvailable, /proc/meminfo) is sampled every 2 s.  The driver prints no
# h_min, so it is derived from its dt line: h_min = dt p^2 c_p / (cfl C(p)), cfl 0.5,
# c_p = 6000 m/s in both decks (the same relation gives 26.81 m for the TPV102 200 m mesh).
for case in "tpv104 tpv104/configs/tpv104_nitsche.toml tpv104/mesh/tpv104_200m.msh" \
            "tpv205 tpv205/configs/tpv205_nitsche.toml tpv205/mesh/tpv2053d_200m.msh"; do
   set -- $case; PROB=$1; DECK=$2; MESH=$3
   echo "=== 200 m load check: $PROB p=2, 144 ranks, tfinal 0.003 s ==="
   OUT=$SCRATCH/nitsche_benchmarks/loadcheck_${PROB}_200m_p2_$SLURM_JOB_ID; mkdir -p $OUT
   ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
   MON=$!
   s=$(date +%s)
   ibrun ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH --order 2 --tfinal 0.003 \
      --output-dir $OUT --log-every 1 > $OUT/driver.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   echo "driver exit $rc in $(( $(date +%s) - s )) s"
   grep -a "tets;" $OUT/driver.log | head -1
   grep -a "  dt = " $OUT/driver.log | head -1
   grep -a "  dt = " $OUT/driver.log | head -1 | awk '{for(i=1;i<=NF;i++){if($i=="dt"&&$(i+1)=="=")dt=$(i+2); if($i=="C(p)")C=$(i+1)}} END{h=dt*4*6000/(0.5*C); printf "implied h_min (min 6V/A) %.2f m; steps to 12 s %d\n", h, int(12/dt+0.999)}'
   grep -a "=== done" $OUT/driver.log || tail -15 $OUT/driver.log
   awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "node memory: peak in use %.1f GB of %.1f GB (MemTotal - MemAvailable, 2 s samples)\n", m, t}' $OUT/meminfo_kB.log
done
sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,Elapsed,MaxRSS,MaxVMSize,NTasks 2>&1 | head -20
echo "(MaxRSS is per task and appears only for srun steps; ibrun with Open MPI may leave it empty, so the /proc/meminfo peak above is the number that decides 144 vs 128 ranks per node)"
echo "total $(( $(date +%s) - T0 )) s"
