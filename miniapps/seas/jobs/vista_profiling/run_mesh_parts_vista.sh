#!/bin/bash
# Shared body of the mesh-setup checks (partition-memory report rev1 options D and E, speed-next S1;
# Chunhui Zhao, 2026-10-22).  Sourced by mp_lts_tpv205_200m_p{1,2}.sbatch, which set PROB, MESH, ORDER
# and VARIANTS, one "name|ranks|deck|flags|tfinal|reference|check" per variant:
#   ranks      all (the job's ranks) or 1 (one rank on the first node: the part file split)
#   flags      the driver flags after --config/--mesh/--order/--tfinal/--output-dir/--log-every (no
#              --lts is added); @BASE is replaced by the job's output folder and @NTASKS by its rank count
#   tfinal     the end time (none for the split, which stops after writing the parts)
#   reference  a gold results folder, @<variant> for an earlier variant of this job, or empty
#   check      bitwise (the job fails unless every station equals the reference's first tfinal s as
#              printed) or roundoff (the differences are printed, nothing is required)
# One allocation runs the variants in turn.  For each it prints the driver's mesh setup line, the
# setup times and the setup memory (VmData per phase, max/mean over ranks), the peak Committed_AS of
# the first node against its CommitLimit (Vista's compute nodes refuse allocations past CommitLimit
# instead of swapping), the stepping time and the station comparison.
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose driver was built at the same
# commit (build_mesh_parts_vista.sbatch, through build_vista.sh, writes miniapps/seas/BUILD_COMMIT_VISTA)
# and whose mesh sha256 is in jobs/vista_profiling/meshes.sha256.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # one thread per MPI rank
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)    # repository root
cd $ROOT
# jobs that start together fetch in the same checkout and can collide on the ref lock ("cannot lock
# ref", job 1075918), so retry after a random pause
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
HEAD=$(git rev-parse HEAD)
cd miniapps/seas
echo "host $(hostname)  nodes $SLURM_JOB_NUM_NODES  ranks $SLURM_NTASKS"
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
[ "$(cat BUILD_COMMIT_VISTA 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_mesh_parts_vista.sbatch"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; build_mesh_parts_vista.sbatch copies it"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" jobs/vista_profiling/meshes.sha256 \
   || { echo "ERROR: $MESH sha256 $SHA is not the gold mesh (jobs/vista_profiling/meshes.sha256)"; exit 2; }
for spec in "${VARIANTS[@]}"; do
   IFS='|' read -r v NR DECK FLAGS TF REF CHECK <<< "$spec"
   case "$REF" in ''|@*) ;; *) [ -d "$REF" ] || { echo "ERROR: gold folder $REF missing"; exit 2; } ;; esac
done
echo "mesh $MESH sha256 $SHA  order $ORDER"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
BASE=$SCRATCH/nitsche_profiling/mp_${PROB}_200m_p${ORDER}_$SLURM_JOB_ID
mkdir -p $BASE
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > $BASE/RUN_INFO.txt
T0=$(date +%s)
rc_all=0
for spec in "${VARIANTS[@]}"; do
   IFS='|' read -r v NR DECK FLAGS TF REF CHECK <<< "$spec"
   FLAGS=${FLAGS//@BASE/$BASE}; FLAGS=${FLAGS//@NTASKS/$SLURM_NTASKS}
   OUT=$BASE/$v; mkdir -p $OUT
   cp "$DECK" $OUT/
   # this (the first) node holds rank 0: its committed memory against the commit limit, and the
   # memory in use
   ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable|^CommitLimit|^Committed_AS/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $OUT/meminfo_kB.log &
   MON=$!
   s=$(date +%s)
   if [ "$NR" = 1 ]; then
      ibrun -n 1 -o 0 ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER \
            --output-dir $OUT $FLAGS > $OUT/driver.log 2>&1
   else
      ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TF \
            --output-dir $OUT --log-every 1000 $FLAGS > $OUT/driver.log 2>&1
   fi
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   [ $rc -eq 0 ] || rc_all=$rc
   echo "=== $v: $NR rank(s), deck $DECK, flags [$FLAGS], tfinal ${TF:-none}: driver exit $rc, $(( $(date +%s) - s )) s in the job"
   grep -a "mesh setup\|mesh parts:\|partition:\|LTS\|setup (wall\|setup memory\|time in ode.Step per rank\|halo\|station output:\|=== done" $OUT/driver.log || tail -20 $OUT/driver.log
   [ $rc -eq 0 ] || tail -20 $OUT/driver.log
   # columns: time, MemTotal, MemAvailable, CommitLimit, Committed_AS (kB, /proc/meminfo order)
   awk '{c=$5/1048576; if(c>m)m=c; l=$4/1048576; u=($2-$3)/1048576; if(u>n)n=u} END{printf "first node: peak Committed_AS %.1f GiB of CommitLimit %.1f GiB; peak in use %.1f GiB (2 s samples)\n", m, l, n}' $OUT/meminfo_kB.log
   [ "$NR" = 1 ] && [ -d "$BASE/parts" ] && echo "parts: $(ls $BASE/parts | wc -l) files, $(du -sh $BASE/parts | cut -f1)"
   if [ -n "$REF" ] && [ $rc -eq 0 ]; then
      case "$REF" in @*) RD=$BASE/${REF#@} ;; *) RD=$REF ;; esac
      echo "--- $v stations against the first $TF s of $RD ($CHECK) ---"
      if [ "$CHECK" = bitwise ]; then
         python3 -I jobs/vista_profiling/compare_stations.py $RD $OUT --tmax $TF --time-tol 2e-9 --groups --bitwise \
            || { rc_all=1; echo "ERROR: $v is not bitwise equal to $RD"; }
      else
         python3 -I jobs/vista_profiling/compare_stations.py $RD $OUT --tmax $TF --time-tol 2e-9 --groups
      fi
   fi
done
du -sh $BASE/*
echo "total $(( $(date +%s) - T0 )) s"
exit $rc_all
