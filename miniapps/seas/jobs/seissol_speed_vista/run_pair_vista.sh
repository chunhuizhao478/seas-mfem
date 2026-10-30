#!/bin/bash
# Shared body of the paired speed-test jobs (SeisSol and MFEM speed test on Vista, plan rev1, Chunhui
# Zhao 2026-10-30) and of the SeisSol layout probe.  Sourced by run_<prob>_p<N>.sbatch and
# probe_layout.sbatch, which set
#   PROB       tpv205 | tpv104 | tpv102
#   MESH       the Gmsh mesh both codes read (path under miniapps/seas)
#   PUML       its pumgen conversion's base name in $WORK/seissol/meshes (tpv2053d_200m, tpv104_200m)
#   ORDER      MFEM polynomial order p; SeisSol runs order o = p + 1 (the July pairing)
#   TFINAL     simulated seconds, the same for both codes
#   MFEM_RUN   1 = run MFEM first (empty in the probe)
#   SS_LAYOUTS SeisSol MPI ranks per node, one run per entry (the probe has three; a run job has the one
#              the probe picked)
# One allocation runs MFEM, then SeisSol, on the same nodes, so the times compare only within a job (S1
# 1102099 against V3 1100111: the same MFEM job 18 % slower on other nodes).  Only stations are written
# by either code (plan section 4).  Prints, per code: the job seconds of each run, MFEM's ode.Step and
# loop lines, SeisSol's "Simulation time" lines and cluster histogram.
# Runs only from a clean pushed checkout (repo CLAUDE.md guard) whose MFEM driver was built at the same
# commit and whose SeisSol binaries and meshes are those of build_speed_test_vista.sbatch.
export LC_ALL=C LANG=C
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)    # repository root
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
HEAD=$(git rev-parse HEAD)
cd miniapps/seas
SEAS=$PWD
SST=$SEAS/jobs/seissol_speed_vista
SEISSOL_VISTA=$WORK/seissol
SO=$(( ORDER + 1 ))
SSBIN=$(find $SEISSOL_VISTA/build-o$SO-elastic-m0-sve128 -maxdepth 1 -type f -perm -u+x -name 'seissol-elastic-o*' ! -name 'proxy*' 2>/dev/null | head -1)
KIT_ENV=$SEISSOL_VISTA/kit/seissol-env-vista.sh
KIT_ENV_SHA256=5df20ce33e9ec4ef9e6deec905ad113a34d5bf934c6021ba2ae7f6f38e5f7209
echo "host $(hostname)  nodes $SLURM_JOB_NUM_NODES  tasks $SLURM_NTASKS  $PROB p$ORDER / o$SO  0 -> $TFINAL s"
# --- inputs: the MFEM build, the mesh, the SeisSol binary, the conversion, the kit environment
if [ -n "$MFEM_RUN" ]; then
   [ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
   [ "$(cat BUILD_COMMIT_VISTA 2>/dev/null)" = "$HEAD" ] \
      || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_VISTA 2>/dev/null)', checkout is $HEAD; rerun build_speed_test_vista.sbatch"; exit 2; }
fi
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" jobs/vista_profiling/meshes.sha256 \
   || { echo "ERROR: $MESH sha256 $SHA is not the gold mesh (jobs/vista_profiling/meshes.sha256)"; exit 2; }
grep -q "^$SHA  $MESH\$" $SEISSOL_VISTA/meshes/MESHES.sha256 \
   || { echo "ERROR: $SEISSOL_VISTA/meshes was not converted from this $MESH"; exit 2; }
PSHA=$(sha256sum $SEISSOL_VISTA/meshes/$PUML.puml.h5 | cut -d' ' -f1)
grep -q "^$PSHA  $PUML.puml.h5\$" $SEISSOL_VISTA/meshes/MESHES.sha256 \
   || { echo "ERROR: $PUML.puml.h5 sha256 $PSHA is not the one build_speed_test_vista.sbatch wrote"; exit 2; }
[ -n "$SSBIN" ] || { echo "ERROR: no SeisSol o$SO binary; run build_speed_test_vista.sbatch"; exit 2; }
grep -q "^   sha256 $(sha256sum $SSBIN | cut -d' ' -f1)\$" $SEISSOL_VISTA/BUILD_MANIFEST-speedtest.txt \
   || { echo "ERROR: $SSBIN is not the binary in BUILD_MANIFEST-speedtest.txt"; exit 2; }
[ "$(sha256sum $KIT_ENV | cut -d' ' -f1)" = $KIT_ENV_SHA256 ] || { echo "ERROR: $KIT_ENV changed"; exit 2; }
echo "mesh $MESH $SHA; $PUML.puml.h5 $PSHA; SeisSol $SSBIN"
# the kit's Python (3.11: tomllib; numpy; PyYAML)
( set +u; source $KIT_ENV > /dev/null 2>&1; $VENV/bin/python3 -I $SST/check_parity.py --problem $PROB --repo $ROOT ) \
   > $SCRATCH/parity_${PROB}_$SLURM_JOB_ID.txt 2>&1 \
   || { cat $SCRATCH/parity_${PROB}_$SLURM_JOB_ID.txt; echo "ERROR: the MFEM and SeisSol decks of $PROB differ (check_parity.py)"; exit 2; }
tail -1 $SCRATCH/parity_${PROB}_$SLURM_JOB_ID.txt
BASE=$SCRATCH/seissol_speed/${TAG:-sst}_${PROB}_200m_p${ORDER}_$SLURM_JOB_ID
mkdir -p $BASE
{ echo "commit $HEAD"; echo "mesh $MESH $SHA"; echo "puml $PUML.puml.h5 $PSHA"; echo "seissol $SSBIN $(sha256sum $SSBIN | cut -d' ' -f1)"
  echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ($SLURM_JOB_NODELIST)"; } > $BASE/RUN_INFO.txt
T0=$(date +%s)
rc_all=0
memlog() { ( while :; do echo "$(date +%s) $(awk '/^MemTotal|^MemAvailable/{print $2}' /proc/meminfo | tr '\n' ' ')"; sleep 2; done ) > $1 & MON=$!; }
memsum() { awk '{u=($2-$3)/1048576; if(u>m)m=u; t=$2/1048576} END{printf "first node memory: peak in use %.1f GB of %.1f GB (2 s samples)\n", m, t}' $1; }

# --- MFEM: central difference, --lts rate2, all cores one rank each (MFEM_RPN ranks per node)
if [ -n "$MFEM_RUN" ]; then
   OUT=$BASE/mfem; mkdir -p $OUT
   DECK=jobs/seissol_speed_vista/mfem_decks/${PROB}_nitsche_sst.toml
   cp $DECK $OUT/
   (
   export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
   module purge; module load TACC
   module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
   module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6
   echo "MFEM modules: $(module -t list 2>&1 | tr '\n' ' ')"
   memlog $OUT/meminfo_kB.log
   s=$(date +%s)
   # 144 ranks per node fill the allocation: ibrun, as every MFEM speed run so far.  Fewer (128 at p = 3, as
   # the p = 3 golds) need an even spread over the nodes, which ibrun -n does not give: mpirun ppr
   if [ $MFEM_RPN = 144 ]; then LAUNCH="ibrun"
   else LAUNCH="mpirun -np $(( SLURM_JOB_NUM_NODES * MFEM_RPN )) --map-by ppr:$MFEM_RPN:node --bind-to core"; fi
   echo "MFEM launcher: $LAUNCH"
   $LAUNCH ./seas_nitsche_dyn_driver --config $DECK --mesh $MESH --order $ORDER \
         --tfinal $TFINAL --output-dir $OUT --log-every 1000 --lts rate2 > $OUT/driver.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   echo "=== MFEM p$ORDER on $(( SLURM_JOB_NUM_NODES * MFEM_RPN )) ranks ($MFEM_RPN per node): driver exit $rc, $(( $(date +%s) - s )) s in the job"
   grep -a "partition:\|LTS\|time in ode.Step per rank\|rank balance\|halo\|setup (wall\|station output:\|loop\|=== done" $OUT/driver.log || tail -20 $OUT/driver.log
   memsum $OUT/meminfo_kB.log
   exit $rc
   ) || rc_all=1
fi

# --- SeisSol: order o = p + 1, ClusteredLTS 2, CFL 0.5, one run per layout (ranks per node x threads)
for RPN in $SS_LAYOUTS; do
   OUT=$BASE/seissol_r$RPN; mkdir -p $OUT/output
   cp $SST/decks/$PROB/* $OUT/
   ln -s $SEISSOL_VISTA/meshes/$PUML.puml.h5 $OUT/$PUML.puml.h5
   # the window: decks hold EndTime = 12.0; a shorter TFINAL (the probe) is set here, in the copy only
   sed -i "s/^\( *EndTime *= *\)12\.0/\1$TFINAL/" $OUT/parameters.par
   grep -q "^ *EndTime *= *$TFINAL" $OUT/parameters.par || { echo "ERROR: EndTime not set to $TFINAL"; rc_all=2; continue; }
   PE=$(( 144 / RPN ))                    # cores per rank; one of them runs SeisSol's communication thread
   (
   set +u
   source $KIT_ENV
   [ "$SEISSOL_ENV_OK" = 1 ] || { echo "ERROR: the kit environment self-check failed"; exit 2; }
   export OMP_NUM_THREADS=$(( PE - 1 )) SEISSOL_COMMTHREAD=1 OMP_PLACES=cores OMP_PROC_BIND=close
   echo "SeisSol modules: $(module -t list 2>&1 | tr '\n' ' ')"
   memlog $OUT/meminfo_kB.log
   s=$(date +%s)
   cd $OUT && mpirun -np $(( SLURM_JOB_NUM_NODES * RPN )) --map-by ppr:$RPN:node:pe=$PE --bind-to core \
        -x OMP_NUM_THREADS -x SEISSOL_COMMTHREAD -x OMP_PLACES -x OMP_PROC_BIND \
        $SSBIN parameters.par > $OUT/seissol.log 2>&1
   rc=$?; kill $MON 2>/dev/null; wait $MON 2>/dev/null
   echo "=== SeisSol o$SO, $RPN ranks per node x ($(( PE - 1 )) + 1) threads, $(( SLURM_JOB_NUM_NODES * RPN )) ranks: exit $rc, $(( $(date +%s) - s )) s in the job"
   grep -a "Simulation time (\|Time spent for initial IO\|Total time spent:\|Minimum timestep\|Theoretical speedup\|Cluster histogram\|^.* : *[0-9]* *, *[0-9]*$\|Load imbalance\|Total calculated HW-FLOP\|error\|Error" $OUT/seissol.log | tail -40
   memsum $OUT/meminfo_kB.log
   exit $rc
   ) || rc_all=1
   [ -n "$SS_SCEC" ] && $SEISSOL_VISTA/venv/bin/python3 -I $SST/to_scec.py --problem $PROB --seissol-dir $OUT/output --out $OUT/scec
done
du -sh $BASE/*
echo "total $(( $(date +%s) - T0 )) s"
exit $rc_all
